#!/usr/bin/perl
# RUN_DSCAN.PL
# WRAPPER CODE FOR DELSCAN
# DATE 10/27/10

use warnings;



my $idir=$ARGV[0];
my $pfile=$ARGV[1];


if ($#ARGV < 1) {usage(); exit;}

check_args(@ARGV);

my @oput=read_control($pfile);
my @files=@{$oput[0]};
my @names=@{$oput[1]};

create_input(\@files,\@names,$idir);

call_dels();

clean_up();




######## FUNCTIONS ############

sub check_args{
    for (my $i=0;$i<2;$i++){
	if (! -e $_[$i]){ die "$_[$i] doesn't exist";}
    }}




sub usage{
print " 

RUN_DSCAN.PL
Usage: input_directory family_file\n\n\n";

}

sub read_control{
my $pfile=$_[0];
my @files;
my @names;

open(PED,"<$pfile") or die "can't find family file $pfile\n";

for (my $i=0;$i<3;$i++){

    $line=<PED>;

    $line =~ s/\t/ /g;

    @tmp=split(" ",$line);

    $files[$i]=$tmp[1];  
    $names[$i]=$tmp[2];
   
}

return (\@files,\@names);

}


sub clean_up{

    for (my $i=1; $i<23;$i++){
	
     if (-e "tmp_chr$i.txt"){
	    system("rm tmp_chr$i.txt");}
    }

    system("rm tmp_pedfile");

}

sub call_dels{

      if (-e "$pfile.dels" ) {  system("rm $pfile.dels");}    
      open(OUT,">$pfile.dels");
      print OUT "chromosome start end nsnp length nme carrier offspring_id\n";
      close (OUT);

 for (my $i=1; $i<23;$i++){
     

     if (-e "tmp_chr$i.txt"  ) {  
      
     system("./delscan tmp_chr$i.txt tmp_pedfile 2 2 $i >> $pfile.dels");
     
     } }


}


sub create_input{
my @files=@{$_[0]};
my @names=@{$_[1]};
my $idir=$_[2];
my $chr;

 
    print STDERR "Making temporary input files: \n";


# MAKE TEMP PEDFILE


   open(OUT,">tmp_pedfile");
   print OUT "1 $names[0] 0 0 1\n";
   print OUT "1 $names[1] 0 0 2\n";
   print OUT "1 $names[2] $names[0] $names[1] 1\n";
   close(OUT);

# MAKE MERGED GT FILES

for ($i=0;$i<=$#files;$i++){ chomp($files[$i]);}


    open(IN1,"<$idir/$files[0]") or die "Can't open input $idir/$files[0]";
    open(IN2,"<$idir/$files[1]");
    open(IN3,"<$idir/$files[2]");
  



    while(<IN1>){
       chomp($line1=$_);
       chomp($line2=<IN2>);
       chomp($line3=<IN3>);
       if ($line1=~/#/){next;}
       
       @d1=split(' ',$line1);
       @d2=split(' ',$line2);
       @d3=split(' ',$line3);


     #Open first ouput file if not already open
       if (! defined $chr ){

       $chr = $d1[1];
       $ofile = "tmp_chr$chr.txt";
       open(OUT,">$ofile") or die "Can't open output $ofile\n";
       print OUT "I id $names[0] $names[0] $names[1] $names[1] $names[2] $names[2]\n";
       print STDERR "Chr $chr\n";
       }

     #Spawn new file if new chromosome
        if ($d1[1] ne $chr){
	   close(OUT);

	   if ($d1[1] eq "Y") { last;}

	   $chr=$d1[1];
	   print STDERR "Chr $chr\n";
           
           $ofile="tmp_chr$chr.txt";
           open(OUT,">$ofile");
           print OUT "I id $names[0] $names[0] $names[1] $names[1] $names[2] $names[2]\n";

       
       }


       if ($chr eq "X") { 


	   $a2[0]=$d2[3];
	   $a2[1]=$d2[3];

	   if (length($d3[3])==2){  @a3=split('',$d3[3]); } else {
	       $a3[0]=$d3[3]; $a3[1]=$d3[3];}
            }
       else {
       @a1=split('',$d1[3]);
       @a2=split('',$d2[3]);
       @a3=split('',$d3[3]);

       }
         
       $oline="M $d1[2] $a1[0] $a1[1] $a2[0] $a2[1] $a3[0] $a3[1]\n";
       $oline=~s/\-/N/g; #Change missing data character from 23andMe format to "N"
       if ($oline =~ /DI/){next;}
       print OUT $oline;
}
close(OUT);


}
