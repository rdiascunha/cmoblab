# Bin data for a histogram
# Modifying Tara Murphy's code
# usage: for data values between 0 and 10, with 5 bins
# gawk -f histogram.awk 0 10 5 data.out data.in  ; range: from 0 to 10 using 5 bins

BEGIN { 
  m = ARGV[1]      # minimum value of x axis
  ARGV[1] = ""
  M = ARGV[2]      # maximum value of x axis
  ARGV[2] = ""
  b = ARGV[3]      # number of bins
  ARGV[3] = ""
  file = ARGV[4]   # output file for plotting
  ARGV[4] = ""

  w = (M-m)/b      # width of each bin
  
  print "number of bins = "b
  print "width of bins  = "w 
  print

  # set up arrays
  for (i = 1 ; i <= b; ++i) {
    n[i] = m+(i*w)        # upper bound of bin
    c[i] = n[i] - (w/2)   # centre of bin
    f[i] = 0              # frequency count
  }
}

{ 
  # bins the data
  for (i = 1; i <= b; ++i)
    if ($1 <= n[i]) {
         ++f[i]
         break        
    } 
}

END {
  # print results to screen
  # and to a file for plotting
  print "bin(centre) = freq"

  for (i = 1; i <= b; ++i) {
    if (f[i] > 0) {
      print "bin("c[i]")", "=", f[i]
      print c[i], f[i] > file
    }
    else {
      print "bin("c[i]")", "=", 0
      print c[i], 0 > file
    }
  }
}
