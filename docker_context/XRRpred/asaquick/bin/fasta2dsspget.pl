#!/usr/bin/perl -w
# Original: Nov. 2008, by Eshel Faraggi (c)
# Modified: 2009,2010,2011,2012,2013
# USE: $0 fasta_file
# Generates a mock dsspget file (only sequence information)

   if ($#ARGV == 0)
     {
     $dsspnm = $ARGV[0];
     open(FL,"$dsspnm") || die "Could not open $!, aborting";
     @ddat = <FL>;
     close(FL);
     }
   else
     {
     @ddat = <STDIN>;
#     die "Use $0 fasta_file to generate a mock dsspget file (asa=-1, angles=-360)\n";
     }
   $i = 0;
   for (@ddat)
     {
     chomp;
     next if (substr($_,0,1) eq ">");
     @tmp = split(//,$_);
     for (@tmp)
       {
       $i++;
       printf("%7d %s %s %s %5d %9.1f %9.1f\n",$i,$_,"-","C",-1,-360,-360);
       }
     }
__END__
