use Math::Complex;

# Hydrophobicity values
my %hyd = (
    R => -4.5,
    K => -3.9,
    N => -3.5,
    D => -3.5,
    Q => -3.5,
    E => -3.5,
    H => -3.2,
    P => -1.6,
    Y => -1.3,
    W => -0.9,
    S => -0.8,
    T => -0.7,
    G => -0.4,
    A => 1.8,
    M => 1.9,
    C => 2.5,
    F => 2.8,
    L => 3.8,
    V => 4.2,
    I => 4.5
);

# Read the UniProt file
my $filename = "protein.swiss";  # Replace with the actual filename
open(my $fh, '<', $filename) or die "Failed to open file: $!";
my $protein_sequence = "";

# Flag to indicate if the SQ line is encountered
my $sequence_started = 0;

# Initialize variables
my @true_labels;
my @predictions;
my @windows;
my @window_positions;
my @start_pos;
my @end_pos;

# Iterate through the file
while (my $line = <$fh>) {
    chomp($line);
    if ($line =~ /^SQ/) {
        # Start of the sequence
        $sequence_started = 1;
    } elsif ($sequence_started && $line !~ /^\/\//) {
        # Inside the sequence
        $line =~ s/\s+//g;  # Remove any whitespace
        # Append the line to the protein sequence variable
        $protein_sequence .= $line;
    }
    if ($line =~ /^FT\s+TRANSMEM\s+(\d+)\.\.(\d+)/) {
        my $start_pos = $1;
        push @start_pos, $start_pos;
        my $end_pos = $2;
        push @end_pos, $end_pos;
    }
    last if $line =~ /^\/\//;  # Stop reading after encountering "//" line
}

# Close the file
close($fh);

# Generate the output CSV file
my $csv_filename = "output.csv";  # Replace with the desired output file name
open(my $csv_fh, '>', $csv_filename) or die "Failed to open CSV file: $!";

print "Enter the window size: ";
my $window_size = <>;
chomp($window_size);

print "Enter the sensitivity: ";
my $sensitivity = <>;
chomp($sensitivity);

# Print the protein sequence
print "Protein sequence: $protein_sequence\n";

my @frame;
for (my $i = int($window_size / 2); $i < length($protein_sequence) - 2; $i += 1) {
    push(@frame, substr($protein_sequence, $i - int($window_size / 2), $window_size));
}

sub calculate_average {
    my @values = @_;
    my $sum = 0;
    foreach my $value (@values) {
        $sum += $value;
    }
    my $average = $sum / scalar(@values);
    $average = sprintf("%.3f", $average);
    return $average;
}

sub hydrophobicity {
    my ($window_size, @array) = @_;
    my @hydrophobicity;
    my $overall_sum = 0;
    my $overall_count = 0;

    for (my $i = 0; $i < scalar(@array); $i++) {
        my $q = $array[$i];
        my @q = split(//, $q);
        for (my $j = 0; $j < scalar(@q); $j++) {
            my $aa = $q[$j];
            if (exists $hyd{$aa}) {
                my $temp = $hyd{$aa};
                push @hydrophobicity, $temp;
                $overall_sum += $hyd{$aa};
                $overall_count++;
            } else {
                push @hydrophobicity, "NA";
            }
        }
    }

    for (my $i = 0; $i <= scalar(@hydrophobicity) - $window_size; $i += $window_size) {
        my @window;
        for (my $j = 0; $j < $window_size; $j++) {
            $window[$j] = @hydrophobicity[$i + $j];
        }
        push @windows, \@window;
    }

    for (my $i = int($window_size / 2); $i <= scalar(@array) + int($window_size / 2); $i++) {
        my $window_start = $i - int($window_size / 2);
        my $window_end = $i + int($window_size / 2);
        my $window_positions = "$window_start-$window_end";
        push @window_positions, $window_positions;
    }

    print $csv_fh "Window Positions   \t   Window Contents   \t   Mean Hydrophobicity\n";

    foreach my $i (0 .. $#windows) {
        my $window_positions = $window_positions[$i];
        my $window = $windows[$i];
        my @amino_acids = $frame[$i];
        my $amino_acids_string = join("", @amino_acids);
        my $mean_hydrophobicity = calculate_average(@$window);
        my $prediction = 0;
        if ($mean_hydrophobicity > 2) {
            $prediction = 1;
        }
        push @predictions, $prediction;
        print $csv_fh "$window_positions                                      \t$amino_acids_string                                     \t$mean_hydrophobicity\n";
    }

    my $overall_average = $overall_sum / $overall_count;
    print $csv_fh "Overall Average Hydrophobicity: $overall_average\n";
}

hydrophobicity($window_size, @frame);

my @window_labels;  # Labels for each window

for (my $i = 0; $i < scalar(@windows); $i++) {
    my $window_positions = $window_positions[$i];
    my ($window_start, $window_end) = split('-', $window_positions);

    # Check if the window overlaps with any known transmembrane region
    my $window_label = 0;  # Default label is 0 (non-transmembrane)

    for (my $j = 0; $j < scalar(@start_pos); $j++) {
        my $region_start = $start_pos[$j];
        my $region_end = $end_pos[$j];
        if ($window_start >= $region_start - $sensitivity && $window_end <= $region_end + $sensitivity) {
            $window_label = 1;  # Window overlaps with transmembrane region
            last;
        }
    }

    push @window_labels, $window_label;
}

# Assign the corresponding labels to @true_labels
@true_labels = @window_labels;

# Calculate true positive (TP), true negative (TN), false positive (FP), and false negative (FN)
my ($TP, $TN, $FP, $FN) = (0, 0, 0, 0);

for (my $i = 0; $i < scalar(@true_labels); $i++) {
    my $true_label = $true_labels[$i];
    my $prediction = $predictions[$i];

    if ($true_label == 1 && $prediction == 1) {
        $TP++;
    } elsif ($true_label == 0 && $prediction == 0) {
        $TN++;
    } elsif ($true_label == 0 && $prediction == 1) {
        $FP++;
    } elsif ($true_label == 1 && $prediction == 0) {
        $FN++;
    }
}

# Calculate accuracy
my $accuracy = ($TP + $TN) / ($TP + $FP + $TN + $FN);

print $csv_fh "True Positives (TP): $TP\n";
print $csv_fh "True Negatives (TN): $TN\n";
print $csv_fh "False Positives (FP): $FP\n";
print $csv_fh "False Negatives (FN): $FN\n";
print $csv_fh "Accuracy: $accuracy\n";
my $MCC = (($TP * $TN) - ($FP * $FN)) / sqrt(($TP + $FP) * ($TP + $FN) * ($TN + $FP) * ($TN + $FN));

print $csv_fh "Matthews Correlation Coefficient (MCC): $MCC\n";

# Close the CSV file
close($csv_fh);

print "Output written to $csv_filename\n";

