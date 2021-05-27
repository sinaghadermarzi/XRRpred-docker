#!/usr/bin/perl -w
# Nov. 2008, by Eshel Faraggi
# USE: /data3/efaraggi/tp.server/bin/tppred/spineX.pl seq_name
# Send error messages to the user, not system log
#   open(STDERR,'<&STDOUT'); $| = 1;   

   if (($#ARGV == 0)||($#ARGV == 1))
     {
     $dsspnm = $ARGV[0];
     open(FL,"$dsspnm") || die "Could not open $!, aborting";
#     @ddat = <FL>;
     @ddat = `dssp $dsspnm`;
     close(FL);
     if ($#ARGV == 1)
       {
       $chain = $ARGV[1];
       }
     else
       {
       $chain = "no";
       }
     }
   else
     {
#     $chain = "no";
#     @ddat = <STDIN>;
     die "Get cryst_index, aa, ss8, ss3, asa, phi, psi values from PDB file\nUse \"$0 pdb_file [chain_letter]\" \n\n";
     }

   if (substr($ddat[0],0,4) eq "====")		# remove extras from beginning of dssp file
     {
     for (@ddat)
       {
       @tmp = split(/[\ \t\n]+/,$ddat[0]);
       shift(@ddat);
       last if ("$tmp[1]" eq "#");
       }
     }

   for (@ddat)			# read chain data
     {
     if ((substr($_,11,1) eq $chain)||($chain eq "no"))
       {
       push(@ddatc, $_);
       }
     }

   for (@ddatc)			# read dssp sequence
     {
     $ss8 = substr($_,16,1);
     @tmp = split;
     if ( ($ss8 eq "G") || ($ss8 eq "H") || ($ss8 eq "I") )
       {
       $ss3 = "H";
       }
     elsif ( ($ss8 eq "B") || ($ss8 eq "E") )
       {
       $ss3 = "E";
       }
     else
       {
       $ss3 = "C";
       }
     $ss8 = "-" if ($ss8 eq ' ');
     $asa = substr($_,35,3);

     @tempval =  split(/[,]/,$_);
     $phi = substr($tempval[4], 24 , 6);
     $psi = substr($tempval[4], 30 , 6);
     $aas = substr($_,13,1);
     if ($aas =~ m/[a-z]/)
       {
       $aas = "C";
       }

     printf "%7s %s %s %s %5d %9.1f %9.1f\n", $tmp[1], $aas, $ss8, $ss3, $asa, $phi, $psi;
     }

__END__

