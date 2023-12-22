# J Ludwig 2016
# My first ever script # Used Perl to the Rescue book to teach myself to code. Some code may be copied from the book
# Used to find microsatellites in fasta files

#!/usr/bin/perl
use strict; use warnings;

print "Type name of file containing ID information\n";
my $IDfile = <>;
chomp $IDfile;
open(my $inID, "< $IDfile") or die $!;
my %usatindex;
my %index;
my @add;
my $idordercount = 0;
my @usat_hold;
while (<$inID>){
	        my $seqID; 
		my $SP;
		($seqID, $SP) = split (/\t/,);
	        $idordercount +=1;
	        my $card = {
		            idcount => "$idordercount",
		            seqid   => "$seqID",
		            sp      => "$SP"
		           };
      	       push @add, $card;
	      }
foreach my $card (@add) {
			 $index{$card->{idcount}} = $card;
			}
#my $thingindex =  $index{5}->{seqid};
#print $thingindex;
#print "index ends \n";
print "Type name of file containing Fasta sequences\n";

my $Fastafile = <>;
chomp $Fastafile;
open(my $inFasta, "< $Fastafile") or die $1;
my $fastaseq; 
my $fastaseq_oneline;
my $fastaordercount = 1;
local $/ = '>';
while (<$inFasta>) {
		    chomp;
		    /\w+.+?\n(.+)/s and $fastaseq = $1 or next;
		    my $test = $fastaseq;
		    $test =~ tr/\r\n//d;
		    $fastaseq_oneline = $test; 
		    my $store = "N";
		    my $repeatcount = 1;
		    my $check_post_match = "N";
                    my $matchcount = 0;
		    my $position_match;
use re 'eval';
		    while ($fastaseq_oneline =~ /(?<!$store)(?<!$check_post_match)([AGTC]{2})\1(?{ $store = $1; })(?{ $check_post_match = substr ($fastaseq_oneline, $- [ 1 ]+2, 2); })(?{ $position_match = $-[ 1 ]+2; })(?{ until($check_post_match ne $store){$repeatcount++; $check_post_match = substr ($fastaseq_oneline, $position_match+=2, 2);} })/g)
			   {
			    my $fastaorder_match_count = join '', ($fastaordercount,$matchcount);
			    my $pos_motiff_repeat = join ' ', ($-[ 1 ], $1, $repeatcount);
#			    print "fmc $fastaorder_match_count\n";
			    @usat_hold = (
					  {motiff => "$1", position => "$-[ 1 ]", repeats => "$repeatcount", matchcount => "$matchcount", fastaorder => "$fastaordercount", fastaorder_match  => "$fastaorder_match_count", pos_motiff_repeat => "$pos_motiff_repeat"}
				         );
#	print @usat_hold;
print " \n";		
	for my $card (@usat_hold) {
$usatindex{$card->{fastaorder_match}} = $card;
}    
			$repeatcount = 1; $matchcount+=1;
			   }
		     $fastaordercount +=1;
		   }		       



close($inID);
close($inFasta);

#foreach my $card (@usat_hold) {
#				print "$card->{motiff}\n";
#}

#print @usat_hold;
#print "end usat_hold\n";

#foreach my $card (@usat_hold) {
#			       $usatindex{$card->{fastaorder_match}} = $card;
#			      }
#my $acard = $usatindex{118}->{pos_motiff_repeat};
#print $acard;
#print "acard above\n";
#print %usatindex;
#print "\n";
#print "\n";


my %count_posmotre2;
my @count_posmotre;

#print "final count bewlo\n";

foreach my $card (%usatindex) {

my $final_count = $usatindex{$card}->{pos_motiff_repeat}; if(defined($final_count)){push @count_posmotre, $final_count}
#print "hi \n";
#my $final_count = $usatindex{$card}->{pos_motiff_repeat};
# 		           my $final_count = $usatindex{$card}->{pos_motiff_repeat}; 
#print "$final_count\n";
#			    push @count_posmotre, $final_count;
#$card++;
#print @count_posmotre;
#print "@count_posmotre[1]\n";
#			    unless (grep (/$final_count/, @count_posmotre)){push @count_posmotre, $final_count}
			      }
#print @count_posmotre;
foreach my $count (@count_posmotre) {
				     $count_posmotre2{$count}++;
				    }
foreach my $book (sort keys %count_posmotre2) {
				     print "$book, $count_posmotre2{$book}";
#				     print "Pos Motiff Numrepeats, $count appeared $count times\n";
#rintf "%-31s %s\n", $book, $count_posmotre2{$book};
#print "\n";
#printf $count_posmotre2{$book};
print "\n";
				    }

#my $thing;
#$thing = $usatindex{1}->{repeats};
#print $thing;
#print "\n";


__END__

my @data_hold;
for(my $cardcount = 1; $cardcount = $idordercount; $cardcount++){
								 foreach (%index) {
		  								   my $id = $index{$cardcount}->{seqid};
										   my $sp = $index{$cardcount}->{sp};
										  }
								 foreach (%usatindex) {
		      								       my $motiff   = $usatindex{$cardcount}->{motiff};
										       my $position = $usatindex{$cardcount}->{position};
										       my $repeats  = $usatindex{$cardcount}->{repeats};
										      }
								 my $card = {
									     seqid    => "$id",
									     sp       => "$sp",
									     motiff   => "$motiff",
									     position => "$position",
									     repeats  => "$repeats"
									    };
								 push @data_hold, $card;
								}
my %data;
foreach my $card (@data_hold) {
			       $data{$card->{seqid}} = $card;
			      }
foreach my $card (%data) {
		 my $sp_check = $data{$card->{sp}};



my $thing = 0;
$thing = $usatindex{10}->{repeats};
print $thing;
print "\n";
