CDF      
      time       lon       lat             CDI       @Climate Data Interface version 2.0.5 (https://mpimet.mpg.de/cdi)   Conventions       CF-1.6     history      �Sun Jun 16 22:18:05 2024: cdo -s -f nc -setcalendar,360_day -yearmean -fldmean -mergetime -apply,-selname,rlut [ ./output/MOST_PLA.0191.nc ./output/MOST_PLA.0192.nc ./output/MOST_PLA.0193.nc ./output/MOST_PLA.0194.nc ./output/MOST_PLA.0195.nc ./output/MOST_PLA.0196.nc ./output/MOST_PLA.0197.nc ./output/MOST_PLA.0198.nc ./output/MOST_PLA.0199.nc ./output/MOST_PLA.0200.nc ./output/MOST_PLA.0201.nc ./output/MOST_PLA.0202.nc ./output/MOST_PLA.0203.nc ./output/MOST_PLA.0204.nc ./output/MOST_PLA.0205.nc ./output/MOST_PLA.0206.nc ./output/MOST_PLA.0207.nc ./output/MOST_PLA.0208.nc ./output/MOST_PLA.0209.nc ./output/MOST_PLA.0210.nc ./output/MOST_PLA.0211.nc ./output/MOST_PLA.0212.nc ./output/MOST_PLA.0213.nc ./output/MOST_PLA.0214.nc ./output/MOST_PLA.0215.nc ./output/MOST_PLA.0216.nc ./output/MOST_PLA.0217.nc ./output/MOST_PLA.0218.nc ./output/MOST_PLA.0219.nc ./output/MOST_PLA.0220.nc ] gmrlut_r3x2_yrmean.nc
Fri May 10 14:02:46 2024: cdo -f nc -s merge spectral56997.nc gaussian56997.nc tempfile56997.nc
Fri May 10 14:02:27 2024: cdo -t ./param56997.tab -f nc -s -sp2gp -dv2uv -setzaxis,zgrid56997.txt -setgrid,t42 grid_56997_02.srv spectral56997.nc     	frequency         year   CDO       @Climate Data Operators version 2.0.5 (https://mpimet.mpg.de/cdo)         time                standard_name         time   units         hours since 191-6-30 00:00:00      calendar      360_day    axis      T               P   lon                standard_name         	longitude      	long_name         	longitude      units         degrees_east   axis      X               @   lat                standard_name         latitude   	long_name         latitude   units         degrees_north      axis      Y               H   rlut                      	long_name         Top Thermal Radiation      units         W/m2   code         �            X                @(      �iN>@��     �i��@��     �i�[@�S     �i��@��    �f� @��    �jM�@�Q�    �i�!@퉀    �iў@���    �i��@���    �i�*@��    �i�'@�4�    �i�.@�P�    �i��@�l�    �iѠ@���    �i�X@���    �i��A �`    �i�rA�`    �i}�A�`    �i�<A
`    �i�A`    �i�RA&`    �i�|A4`    �i��AB`    �i��A	P`    �i�5A
^`    �i�MAl`    �iV�Az`    �i�JA�`    �i�vA�`    �i��