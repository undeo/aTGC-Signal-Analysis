
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.0 --par1_maximum 0.3 --expected --toys 50
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.0 --par1_maximum -0.3 --expected --toys 50
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.0 --par1_maximum 0.05 --expected --toys 50
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.0 --par1_maximum -0.1 --expected --toys 50
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.0 --expected --toys 50
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.1 --phi_value 45. --expected --toys 50
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.1 --phi_value 120. --expected --toys 50
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.1 --phi_value 200. --expected --toys 50
# python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum -0.3 --par2_minimum 0.0 --par2_maximum -0.1 --phi_value 200. --expected --toys 50


import os
#os.system("date")

import sys
import math

import sys, getopt

def main(argv):
    ws_name = ''
    #   outputfile = ''
    CL = ''
    prec = 0.001
    par = ''
    no_par = ''
    min_val1 = ''
    max_val1 = ''
    min_val2 = ''
    max_val2 = ''
    phi_par = 0.
    exp = False
    noUnc = False
    toys=500.
    
    try:
        opts, args = getopt.getopt(argv,"hi:C:p:a:x:X:y:Y:j:t:en",["inputfile=","CLvalue=","precision=","parameter=","par1_minimum=","par1_maximum=","par2_minimum=","par2_maximum=","phi_value=","toys=",'expected','noUnc'])
    except getopt.GetoptError:
        print 'Error'
        print 'python run_python_1or2D.py --inputfile <input ws name> --CLvalue <CL value for limit> --precision <precision value> --parameter <lZ || dkg || 2D> --par1_minimum <par1 min value> --par1_maximum <par1 max value> --par2_minimum <par2 min value> --par2_maximum <par2 max value> --phi_value <phi value in deg, angle of limit line> --expected --toys <number of toys>'
        print 'Examples:'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.0 --par1_maximum 0.3 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.0 --par1_maximum -0.3 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.0 --par1_maximum 0.1 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.0 --par1_maximum -0.1 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.0 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.1 --phi_value 45. --expected --toys 50'

        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':

            print 'Examples:'
            print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.0 --par1_maximum 0.3 --expected --toys 50'
            print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.0 --par1_maximum -0.3 --expected --toys 50'
            print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.0 --par1_maximum 0.1 --expected --toys 50'
            print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.0 --par1_maximum -0.1 --expected --toys 50'
            print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.0 --expected --toys 50'
            print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.1 --phi_value 45. --expected --toys 50'
            sys.exit()
        elif opt in ("-i", "--inputfile"):
            ws_name = arg
        elif opt in ("-C", "--CLvalue"):
            CL = float(arg)
        elif opt in ("-p", "--precision"):
            prec = float(arg)
        elif opt in ("-a", "--parameter"):
            par = arg
        elif opt in ("-x", "--par1_minimum"):
            min_val1 = float(arg)
        elif opt in ("-X", "--par1_maximum"):
            max_val1 = float(arg)
        elif opt in ("-y", "--par2_minimum"):
            min_val2 = float(arg)
        elif opt in ("-Y", "--par2_maximum"):
            max_val2 = float(arg)
        elif opt in ("-t", "--toys"):
            toys = int(arg)
        elif opt in ("-j", "--phi_value"):
            phi_par = float(arg)
        elif opt in ("-e", "--expected"):
            exp = True
        elif opt in ("-n", "--noUnc"):
            noUnc = True
         
            
    
    if ws_name=='' or CL=='' or par==''  or min_val1==''  or max_val1=='':
        
        print 'Error -> did not define all necesary parameters (--inputfile, --CLvalue, --parameter, --par1_minimum, --par1_maximum) '
        print 'Error'
        print 'python run_python_1or2D.py --inputfile <input ws name> --CLvalue <CL value for limit> --precision <precision value> --parameter <lZ || dkg || 2D> --par1_minimum <par1 min value> --par1_maximum <par1 max value> --par2_minimum <par2 min value> --par2_maximum <par2 max value> --phi_value <phi value in deg, angle of limit line> --expected --toys <number of toys>'
        print 'Examples:'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.0 --par1_maximum 0.3 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.0 --par1_maximum -0.3 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.0 --par1_maximum 0.1 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.0 --par1_maximum -0.1 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.0 --expected --toys 50'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.0 --par1_maximum 0.3 --par2_minimum 0.0 --par2_maximum 0.1 --phi_value 45. --expected --toys 50'
        
        sys.exit(2)


    if par=='lZ':
        no_par='dkg'
    if par=='dkg':
        no_par='lZ'
    if no_par=='' and not par=='2D':
        print "ERROR: did not recognize parameter name!"
        sys.exit()


    if par=="2D" and (min_val2=='' or max_val2==''):
        print 'Error -> did not define all necesary parameters for 2D (--inputfile, --CLvalue, --precision, --parameter, --par1_minimum, --par1_maximum, --par2_minimum, --par2_maximum) '
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter lZ --par1_minimum 0.05 --par1_maximum 0.3'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter dkg --par1_minimum 0.05 --par1_maximum 0.1'
        print 'python run_python_1or2D.py -i Example_ZZ_SMasBKG.root -C 0.95 -p 0.01 -a 2D -x 0.05 -X 0.3 -y 0.01 -Y 0.1 -j 45.'
        print 'python run_python_1or2D.py --inputfile Example_ZZ_SMasBKG.root --CLvalue 0.95 --precision 0.01 --parameter 2D --par1_minimum 0.05 --par1_maximum 0.3 --par2_minimum 0.01 --par2_maximum 0.1 --phi_value 45.'       
        sys.exit(2)
       
    
    print "-->  running FC limits for: ws_name= ", ws_name, " CL= ",CL," precision= ", prec," par= ", par , "  min_par= ", min_val1, "  max_par= ", max_val1, " toys= ", toys
    if par=="2D":
        print " --> 2D limits,  min_par2= ", min_val2, "  max_par2= ", max_val2
    
    # deg->rad
    phi_par=2*math.pi*phi_par/360.
    min_val1=min_val1*math.cos(phi_par)
    max_val1=max_val1*math.cos(phi_par)
    
    if (par=="2D"):
        min_val2=min_val2*math.sin(phi_par)
        max_val2=max_val2*math.sin(phi_par)
        print "    running 2D limit ", "  min_par2= ", min_val2, "  max_par2= ", max_val2, "  phi_par= ", phi_par
    else:
        min_val2=0.
        max_val2=0.
        phi_par=0.
        
    
    file_output = open('FC_results_%s_CL%s_prec%s_par%s_par1min%s_par1max%s_par2min%s_par2max%s_phi%s.txt'%(ws_name,CL,prec,par,min_val1,max_val1,min_val2,max_val2,phi_par), 'w')
   
    file_output.write('/////////////////////////////////////////////////////////////////////////\n"')
    if (par=="2D"):
        file_output.write('//  running FC limits for: ws_name=%s   CL= %s    precision= %s    par= %s   min_par1= %s   max_par1=%s   min_par2= %s   max_par2=%s  phi=%s\n"'%(ws_name,CL,prec,par,min_val1,max_val1,min_val2,max_val2,phi_par))
    else:
        file_output.write('//  running FC limits for: ws_name=%s   CL= %s    precision= %s    par= %s   min_par= %s   max_par=%s\n"'%(ws_name,CL,prec,par,min_val1,max_val1))
         

    # read CLsb and CLsb_err from file containing combine output
    def returnCLs(fname):
    
        CLsb_tmp=-1000
        CLsberr_tmp=-1000
        datafile = file(fname)
        combine_output_ok=False
        for line in datafile:
            if "Hybrid New" in line:
                combine_output_ok=True
#                print "found HN in line: ", line 
            if combine_output_ok:
                if "CLsplusb" in line:
                    index=line.index("+/-")
                    #                print "index= ",index
                    if not index>0:
                        print " CLsplusb not found in output file! -> EXIT"
                        sys.exit();
                        
                    print "cls value ", line
                    data = line.split()
                    CLsb_tmp=float(data[2])
                    CLsberr_tmp=float(data[4])
                    return (CLsb_tmp,CLsberr_tmp)

    
    combine_add=''
    name_add=''
    if exp:
        combine_add+=' -t -1 --expectSignal=1'
        name_add+="EXP"
    if noUnc:
        combine_add+=' -S 0'
        name_add+="noUnc"

    print "\n runing combine for min and max points ..."
    #running combine for max and min parameter values
    if not par=="2D":
        print "  running 1D limits ... "
        file_output.write('1D limits\n')
        file_output.write('par(%s) \t CLsb \t CLsb_err \t output_file\n'%(par))
        
        print "    combine -n _%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s --redefineSignalPOIs %s --freezeNuisances dg1,r,%s %s > out_%s_%s_test_FC_%s_min"%(par,min_val1,name_add,ws_name,int(toys),par,min_val1,par,no_par,combine_add,par,min_val1,name_add)
        os.system("combine -n _%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s --redefineSignalPOIs %s --freezeNuisances dg1,r,%s %s > out_%s_%s_test_FC_%s_min"%(par,min_val1,name_add,ws_name,int(toys),par,min_val1,par,no_par,combine_add,par,min_val1,name_add))
        
        print "    combine -n _%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s --redefineSignalPOIs %s --freezeNuisances dg1,r,%s %s> out_%s_%s_test_FC_%s_max"%(par,max_val1,name_add,ws_name,int(toys),par,max_val1,par,no_par,combine_add,par,max_val1,name_add)
        os.system("combine -n _%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s --redefineSignalPOIs %s --freezeNuisances dg1,r,%s %s > out_%s_%s_test_FC_%s_max"%(par,max_val1,name_add,ws_name,int(toys),par,max_val1,par,no_par,combine_add,par,max_val1,name_add))
        
        CLsb_min, CLsberr_min =returnCLs("out_%s_%s_test_FC_%s_min"%(par,min_val1,name_add))
#        print " -> output min : ",CLsb_min," ",CLsberr_min
        
        CLsb_max, CLsberr_max =returnCLs("out_%s_%s_test_FC_%s_max"%(par,max_val1,name_add))
#        print " -> output max : ",CLsb_max," ",CLsberr_max
        
        file_output.write('%s \t %s \t %s \t higgsCombine_%s_%s_%s_.HybridNew.mH120.root\n'%(min_val1,CLsb_min,CLsberr_min,par,min_val1,name_add))
        file_output.write('%s \t %s \t %s \t higgsCombine_%s_%s_%s_.HybridNew.mH120.root\n'%(max_val1,CLsb_max,CLsberr_max,par,max_val1,name_add))
        
        
    else:
        print "  running 2D limits ..."
        file_output.write('2D limits\n')
        file_output.write('lZ \t dkg \t CLsb \t CLsb_err \t output_file\n')
        
        print "    combine -n _%s_%s_%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s,%s=%s --redefineSignalPOIs %s,%s --freezeNuisances dg1,r %s > out_%s_%s_%s_%s_%s_%s_test_FC_%s_min"%("lZ",min_val1,"dkg",min_val2,name_add,ws_name,int(toys),"lZ",min_val1,"dkg",min_val2,"lZ","dkg",combine_add,"lZ",min_val1,"dkg",min_val2,"phi",phi_par,name_add)
        os.system("combine -n _%s_%s_%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s,%s=%s --redefineSignalPOIs %s,%s --freezeNuisances dg1,r %s > out_%s_%s_%s_%s_%s_%s_test_FC_%s_min"%("lZ",min_val1,"dkg",min_val2,name_add,ws_name,int(toys),"lZ",min_val1,"dkg",min_val2,"lZ","dkg",combine_add,"lZ",min_val1,"dkg",min_val2,"phi",phi_par,name_add))
    
        print "    combine -n _%s_%s_%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s,%s=%s --redefineSignalPOIs %s,%s --freezeNuisances dg1,r %s > out_%s_%s_%s_%s_%s_%s_test_FC_%s_max"%("lZ",max_val1,"dkg",max_val2,name_add,ws_name,int(toys),"lZ",max_val1,"dkg",max_val2,"lZ","dkg",combine_add,"lZ",max_val1,"dkg",max_val2,"phi",phi_par,name_add)
        os.system("combine -n _%s_%s_%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s,%s=%s --redefineSignalPOIs %s,%s --freezeNuisances dg1,r %s > out_%s_%s_%s_%s_%s_%s_test_FC_%s_max"%("lZ",max_val1,"dkg",max_val2,name_add,ws_name,int(toys),"lZ",max_val1,"dkg",max_val2,"lZ","dkg",combine_add,"lZ",max_val1,"dkg",max_val2,"phi",phi_par,name_add))
        
        CLsb_min, CLsberr_min =returnCLs("out_%s_%s_%s_%s_%s_%s_test_FC_%s_min"%("lZ",min_val1,"dkg",min_val2,"phi",phi_par,name_add))
        print " -> output min : ",CLsb_min," ",CLsberr_min
        
        CLsb_max, CLsberr_max =returnCLs("out_%s_%s_%s_%s_%s_%s_test_FC_%s_max"%("lZ",max_val1,"dkg",max_val2,"phi",phi_par,name_add))
        print " -> output max : ",CLsb_max," ",CLsberr_max
        
        file_output.write('%s \t %s \t %s \t %s \t higgsCombine_%s_%s_%s_%s_%s_.HybridNew.mH120.root\n'%(min_val1,min_val2,CLsb_min,CLsberr_min,"lZ",min_val1,"dkg",min_val2,name_add))
        file_output.write('%s \t %s \t %s \t %s \t higgsCombine_%s_%s_%s_%s_%s_.HybridNew.mH120.root\n'%(max_val1,max_val2,CLsb_max,CLsberr_max,"lZ",max_val1,"dkg",max_val2,name_add))


    # make sure that CL_min<CL<CL_max
    print 'CLsb_max= ',CLsb_max,' CLsb_min= ',CLsb_min,'  CL= ',CL
    if not ((CLsb_max>(1-CL) and (1-CL)>CLsb_min) or (CLsb_max<(1-CL) and (1-CL)<CLsb_min) ):
        print 'Error: Desired CLsb (%s) value is not between CLsb values of min (%s) and max (%s) points'%(1-CL,CLsb_min,CLsb_max)
        sys.exit()


    diff=100.

    # run for point in the middle, read CLsb value and decide which way to go for next point
    # run until found limit within precision
    while ( math.fabs(diff) > prec ):
        middle_val1=(max_val1+min_val1)/2.
        middle_val2=(max_val2+min_val2)/2.
        #     print "middle_val1= ", middle_val1
        print "\n\n  --->> running middle point: ", middle_val1," ",middle_val2
        
        if not par=="2D":
            print " running 1D limits"
            
            print "    combine -n _%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s --redefineSignalPOIs %s --freezeNuisances dg1,r,%s %s > out_%s_%s_test_FC_%s"%(par,middle_val1,name_add,ws_name,int(toys),par,middle_val1,par,no_par,combine_add,par, middle_val1,name_add)
            os.system("combine -n _%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s --redefineSignalPOIs %s --freezeNuisances dg1,r,%s %s > out_%s_%s_test_FC_%s"%(par,middle_val1,name_add,ws_name,int(toys),par,middle_val1,par,no_par,combine_add,par, middle_val1,name_add))
            CLsb_middle, CLsberr_middle =returnCLs("out_%s_%s_test_FC_%s"%(par,middle_val1,name_add))
            file_output.write('%s \t %s \t %s \t higgsCombine_%s_%s_%s_.HybridNew.mH120.root\n'%(middle_val1,CLsb_middle,CLsberr_middle,par,middle_val1,name_add))
            
        else:
            print " running 2D limits"
            print "    combine -n _%s_%s_%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s,%s=%s --redefineSignalPOIs %s,%s --freezeNuisances dg1,r %s > out_%s_%s_%s_%s_%s_%s_test_FC_%s"%("lZ",middle_val1,"dkg",middle_val2,name_add,ws_name,int(toys), "lZ" ,middle_val1,"dkg" ,middle_val2, "lZ","dkg",combine_add, "lZ",middle_val1,"dkg" ,middle_val2,"phi",phi_par,name_add)
            os.system("combine -n _%s_%s_%s_%s_%s_ %s -M HybridNew --freq --testStat=PL --rule=CLsplusb --toysH=%s --singlePoint %s=%s,%s=%s --redefineSignalPOIs %s,%s --freezeNuisances dg1,r %s > out_%s_%s_%s_%s_%s_%s_test_FC_%s"%("lZ",middle_val1,"dkg",middle_val2,name_add,ws_name,int(toys), "lZ" ,middle_val1,"dkg" ,middle_val2, "lZ","dkg",combine_add, "lZ",middle_val1,"dkg" ,middle_val2,"phi",phi_par,name_add))
            CLsb_middle, CLsberr_middle =returnCLs("out_%s_%s_%s_%s_%s_%s_test_FC_%s"%( "lZ",middle_val1,"dkg",middle_val2,"phi",phi_par,name_add))
            file_output.write('%s \t %s \t %s \t %s \t higgsCombine_%s_%s_%s_%s_%s_.HybridNew.mH120.root\n'%(middle_val1,middle_val2,CLsb_middle,CLsberr_middle,"lZ",middle_val1,"dkg",middle_val2,name_add))
            
            
            #    CLsb_middle, CLsberr_middle =returnCLs("out_1D_test_FC_%s_%s"%(middle_val1))
        print " -> output: ",CLsb_middle," ",CLsberr_middle
        if (CLsb_middle > (1-CL)):
            print " -> CLsb_middle>%s  --> setting min val1ue"%(1-CL)
            min_val1=middle_val1
            min_val2=middle_val2
        elif (CLsb_middle < (1-CL)):
            print " -> CLsb_middle<%s  --> setting max val1ue"%(1-CL)
            max_val1=middle_val1
            max_val2=middle_val2
        else:
            print " -> CLsb_middle=%s  --> found the limit"%(1-CL)
            CLsb_middle=1-CL

        diff=CLsb_middle + CL - 1
        
    if not par=="2D":
        print "DONE: limit ", CL ," CL at  %s "%(par), middle_val1 ,"   -> CLsb= ",CLsb_middle
    else:
        print "DONE: limit ", CL ," CL at lZ= ", middle_val1,"  dkg= ",middle_val2 ,"   -> CLsb= ",CLsb_middle
    file_output.write("done, limit %s CL at %s %s   -> CLsb= %s "%(CL, middle_val1,middle_val2,CLsb_middle))

if __name__ == "__main__":
   main(sys.argv[1:])


