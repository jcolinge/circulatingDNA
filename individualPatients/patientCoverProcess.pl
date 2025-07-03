#!/usr/bin/perl

use strict;
use Getopt::Long;
use Errno qw(EAGAIN);
use POSIX ":sys_wait_h";
use File::Basename;

my $baseExec = '/data2/jcolinge/fragmentomics/atlas';
my $chrSizeFile = "$baseExec/pos_chromosomes.tsv";
my $outFolder = 'indiv-profiles';

my $nRLCover = 100;
my $atlasFile = "cris-healthy-wps-peaks/compiled-selection.txt";
my $centroFile = "$baseExec/centromere-region.txt";

my $nproc = 1;

my $help;
my $force;
my $noDeleteTmp;
my $testMode;
my $synchro;
my $withStartEnd;
my $chosenChr;

if (!GetOptions('help' => \$help,
                'h' => \$help,
                'testmode' => \$testMode,
                'force' => \$force,
                'synchro' => \$synchro,
                'nodeltmp' => \$noDeleteTmp,
                'nproc=i' => \$nproc,
                'outfolder=s' => \$outFolder,
                'atlas=s' => \$atlasFile,
                'chrsizes=s' => \$chrSizeFile,
                'centro=s' => \$centroFile,
                'exec=s' => \$baseExec
               ) || defined($help) || ($nproc<1)){
  print STDERR "patientCoverProcess.pl [options] bed_files
  
  Extract individual coverage information on the basis of an atlas of nucleosome positions
  
  General options:
    --help
    -h
    -testmode            generates the log file without actually executing the commands
    -force               force processing even if output files exist already
    -nodeltmp            does not delete files created at intermediary steps
    -synchro             all patient processings can start at the same time, default is every 3 min to limit file access conflicts
    --nproc=i            number of processors (default $nproc)
    --outfolder=s	     provide an alternative output folder name (default '$outFolder')
    --exec=s             folder in which are located the executables (default '$baseExec')
    --atlas=s            file containing the atlas of chromosomal positions (default $atlasFile)
    --centro=s           provide an alternative file to define centromer positions (default '$centroFile')
    --chrsizes=s         provide an alternative file to define chromosome sizes (default '$chrSizeFile')
 
  \n";
  exit(0);
}

my $coverageExec = "$baseExec/coverage";
my $smoothExec = "$baseExec/sgfilter";
my $readCoverExec = "$baseExec/readCoverage";

my $failedCommands;

# prepare clean exit in case of interruption
$SIG{INT} = \&signalHandler;
$SIG{TERM} = \&signalHandler;
$SIG{CHLD} = sub{wait};
use Proc::ProcessTable;
sub signalHandler{
  # kill all child processes and exits reporting a message
  my $parent = $$;
  my $procTable = Proc::ProcessTable->new();
  foreach my $proc (@{$procTable->table()}){
    kill(15,$proc->pid) if ($proc->ppid == $parent);
  }
  die("Interrupted by a signal: $!\n");
}


# read centromere coordinates (to get chromosome names!!)
printLog("Loading centromere\t$centroFile");
open(C,$centroFile) || die("Failed loading centromere: $centroFile");
my @centro;
$_ = <C>; # skip column headers
while(<C>){
  next if(/^Y/);
  s/[\n\r]$//;
  my @line = split(/\s+/);
  push(@centro,[$line[0],$line[1],$line[2]]); # chr start stop
}
close(C);


# compile the list of all the bed files to process
my (@bed,@stdin);
if ($ARGV[0] eq '-'){
  # Read from STDIN
  @stdin = <STDIN>;
  shift;
}
foreach (@ARGV, @stdin){
  foreach my $file (glob($_)){
    if ((-f $file) && (index($file,'/') == 0)){
      push(@bed,$file);
    }
  }
}


# launches the fragment processing programs ====================================

my $pid;
my $perproc = (@bed % $nproc == 0) ? @bed/$nproc : int(@bed/$nproc)+1;
for (my $fork = 0; $fork < $nproc; $fork++){
  FORK:{
    if ($pid = fork()){
      # parent
    }
    elsif (defined($pid)){
      # child
      for (my $i = $fork*$perproc; ($i < ($fork+1)*$perproc) && ($i < @bed); $i++){
        my $bed = basename($bed[$i]);
        my $sample = substr($bed,0,index($bed,'.'));
        my @dirs = split(/\//,dirname($bed[$i]));
        my $dir = join('/',@dirs[0 .. ($#dirs)]);

        # creates the sample folder
        my $folder = "$outFolder/$sample"; 
        unless (-d $folder){
          my $cmd = "mkdir $folder";
          executeCommand($cmd);
        }
        
        # launches the programs
        if ($force || !(-f "$folder/chrX_after_selection_cover.txt")){
          sleep($fork*180) unless(defined($synchro));
          
          # coverage (no starts, ends)
          my $cmd = "$coverageExec -sizerange chroma -chrsizes $chrSizeFile -infolder $dir/ -fileregexp $sample -outfolder $folder/";
          executeCommand($cmd);
          
          # smooth raw coverage
          for (my $k = 0; $k < 23; $k++){
            my $chr = $centro[$k]->[0];
            $cmd = "$smoothExec -sx 0.0 -i $folder/chr$chr"."_coverage.tsv -o $folder/smooth-chr$chr"."_cover_no_norm.txt -nl $nRLCover -nr $nRLCover -m 4";
            executeCommand($cmd);
          }
          
          # extract coverage at atlas locations
          $cmd = "$readCoverExec -atlas $atlasFile -coverfolder $folder/ -suffix _cover_no_norm.txt > $folder/$sample-peak-coverage.txt";
          executeCommand($cmd);

          # clean up
          unless(defined($noDeleteTmp)){
            # remove intermediary files
            $cmd = "rm $folder/*.tsv";
            executeCommand($cmd);
            $cmd = "rm $folder/smooth-*.txt";
            executeCommand($cmd);
          }
        }
      }
      exit(0);
    }
    elsif ($! == EAGAIN){
      # recoverable fork error
      sleep(5);
      redo FORK;
    }
    else{
      printLog("Cannot fork ($fork)\t$!");
      exit(1);
    }
  }
}
# waits for all child processes
my $kid;
do{
  $kid = waitpid(-1,&WNOHANG);
} until ($kid == -1);

printLog("Finished");
print STDERR "Failed nonfatal commands:\n$failedCommands\n" if ($failedCommands);



########################

sub printLog{
  my ($text) = @_;
  my @t = localtime(time());
  printf STDERR "%04d-%02d-%02d\t%02d:%02d:%02d\t$text\n",1900+$t[5],$t[4],$t[3],$t[2],$t[1],$t[0];
} # printLog


sub executeCommand{
  my ($cmd, $fatal) = @_;
  printLog("executing\t$cmd");
  unless ($testMode){
    if (system($cmd)){
      printLog("failed\t$cmd");
      $failedCommands .= "$cmd\n";
      exit($fatal) if (defined($fatal) && ($fatal>0));
    }
    else{
      printLog("done");
    }
  }
} # executeCommand
