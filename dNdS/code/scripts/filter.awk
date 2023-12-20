NR > 5 {
  if ($6 < 0.01 || $6 > 2) next;
  # Accumulate sum and count of valid values
  sum += $6;
  count++;
}
END {
  # Calculate and print average if any valid values exist
  if (count > 0) {
    average = sum / count;
    printf(average);
  } else {
    print "NA"
  }
}
