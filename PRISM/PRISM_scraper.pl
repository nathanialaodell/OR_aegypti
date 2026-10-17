#!/usr/bin/perl
use strict;
use warnings;
use DateTime;
my $day;
my $clim_var = var; # MUST CHANGE TO ANY OF `tmax`, `tmin`, `tmean`, etc. (see https://prism.oregonstate.edu/documents/PRISM_downloads_web_service.pdf)
my $base_url = 'https://services.nacse.org/prism/data/get/us/800m';
my $start = DateTime->new( day => 20, month => 8, year => 2025 );
my $stop = DateTime->new( day => 1, month => 10, year => 2025 );
my $format = '?format=bil'; # have to specify that we want bil so we can raster stack in R
while($start <= $stop) {
    $day = $start->strftime('%Y%m%d'); #place date in proper format
    system("curl -L -OJ \"$base_url/$clim_var/$day$format\""); # no slash between day and format for url
    sleep 2; #DO NOT ADJUST
    $start->add(days => 1);
}