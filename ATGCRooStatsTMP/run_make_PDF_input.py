import os
from ROOT import *


print "python make_PDF_input.py --POI cwww"
os.system("python make_PDF_input.py --POI cwww")
print "python make_PDF_input.py --POI ccw"
os.system("python make_PDF_input.py --POI ccw")
print "python make_PDF_input.py --POI cb"
os.system("python make_PDF_input.py --POI cb")
print "python make_PDF_input.py --POI cwww,ccw"
os.system("python make_PDF_input.py --POI cwww,ccw")
print "python make_PDF_input.py --POI ccw,cb"
os.system("python make_PDF_input.py --POI ccw,cb")
print "python make_PDF_input.py --POI cwww,ccw,cb"
os.system("python make_PDF_input.py --POI cwww,ccw,cb")
