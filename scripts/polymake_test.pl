use warnings;
use strict;
use application "polytope";

chomp(my $n_tests = <>);

for (my $t = 0; $t < $n_tests; ++$t) {

    # read the dimensions of the current testcase
    chomp(my $line = <>);
    my @linearray = split(" ", $line);
    my $n = $linearray[0];
    my $d = $linearray[1];

    # read the current testcase
    my @generators;
    for (my $i = 0; $i < $n; ++$i) {
        chomp($line = <>);
        @linearray = split(" ", $line);
        my $target = new Vector<Rational>(@linearray);
        my $source = -$target;
        my $endpoints = ones_vector(2) | (new Matrix<Rational>([$source, $target]));
        my $segment = new Polytope<Rational>(VERTICES=>$endpoints);
        push(@generators, $segment);
    }

    # perform the test
    my $z = minkowski_sum_fukuda(@generators);
    my $n_vertices = $z->N_VERTICES;
    print "$n $d $n_vertices\n";

    # there is always an empty line
    chomp($line = <>);
}
