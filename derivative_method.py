from scipy import optimize
from scipy import interpolate
import matplotlib.pyplot as plt
import pandas as pd
import math
import json

           # This program is designed for the physical chemistry experiment: measuring rate constant of ethyl acetate
           # saponification in basic solutions. In this experiment, conductivity is used to measure concentration and
           # [NaOH] and [EtOAc] is designed to be the same.
           
           # We develop this program based on the derivative method, which set up a linear relationship between derivative of
           # conductivity to time and a polynomial only related to conductivity, so error introduced by determining reaction time
           # can be theoretically eliminated.
           
           # You can use this program to calculate activation energy like this:
           # import derivative_method
           # derivative_method.get_activation_energy(filenames,temperatures,k_0,k_inf,c_0,reportname[optional])
           # parameters: filenames: array, the excel files (conductivity-time raw data) that you need to analyze
           #             temperatures: the temperature that the reaction proceeds, sequence corresponds to filenames
           #             k_0: initial conductivity, measured using NaOH whose concentration is c_0, sequence corresponds to filenames
           #             k_inf: end conductivity, measured using NaOAc whose concentration is c_0, sequence corresponds to filenames
           #             c_0: initial concentration of substrates, sequence corresponds to filenames
           #             reportname[optional]: file name of the analysis result (.json and .png), default is 'report'               


def __init__():

    global sheetname,cols,calibrating_remove_item_count,R_square_threshold,bad_data_threshold,bad_data_initial_value,bad_data_popping_times,molar_gas_constant
    
    print('loading config file')
    
    with open('config.json', 'r', encoding='utf-8') as file:
        config = json.load(file)
    
    sheetname,cols,calibrating_remove_item_count,R_square_threshold,bad_data_threshold,bad_data_initial_value,bad_data_popping_times,molar_gas_constant=\
    config["sheetname"],config["cols"],config["calibrating_remove_item_count"],config["R_square_threshold"],\
    config["bad_data_threshold"],config["bad_data_initial_value"],config["bad_data_popping_times"],config["molar_gas_constant"]



def __conduct_time_func__(x,a,b,c):
    return (a+b*x)/(1+c*x)

def __linear__(x,a,b):
    return a*x+b

def __io__(filepath):
    
    print('loading file '+ filepath + '...')
    
    data=pd.read_excel(filepath,sheet_name=sheetname,usecols=cols)
    data=data.values.tolist()
    time,k_t=[],[]
    for item in data:
        time.append(item[0])
        k_t.append(item[1])
    
   # print(time,k_t)
   # plt.plot(time,k_t)
    
    return time,k_t

    
def get_rate_constant(filepath,k_0,k_inf,c_0):

    __init__()

    a,b=__io__(filepath)
    
    print('Analysing data from '+ filepath + '...')
    
    time,k_t=__data_pretreatment__(a,b)
    
   # plt.plot(time,k_t)
    
    R_2_k_t=0
    
    while R_2_k_t<R_square_threshold:
       __remove_first_items__(k_t,calibrating_remove_item_count)
       __remove_first_items__(time,calibrating_remove_item_count)
       k_t_pred=[]
       pred_args,pcov_1=optimize.curve_fit(__conduct_time_func__,time,k_t)
       for t in time:
           k_t_pred.append(float(__conduct_time_func__(t,pred_args[0],pred_args[1],pred_args[2])))
       R_2_k_t=__determination_coefficient__(k_t,k_t_pred)
       
    #   print(pred_args,k_t,k_t_pred,R_2_k_t)
    
    #plt.plot(time,k_t_pred)
    #plt.plot(time,k_t)
    #plt.show()
    
    k_t_func=interpolate.UnivariateSpline(time,k_t_pred)
    k_t_derivative_func=k_t_func.derivative()
    k_t_derivative=[]
    for t in time:
        k_t_derivative.append(float(k_t_derivative_func(t)))
    
   # plt.plot(time,k_t_derivative)
   # plt.show()
   
    line_x_axis=[]
    for kt in k_t_pred:
        line_x_axis.append((k_inf-kt)**2)
        
   # plt.scatter(line_x_axis,k_t_derivative)
   # plt.show()
    
    k_linear_args,pcov_2=optimize.curve_fit(__linear__,line_x_axis,k_t_derivative)
    
    k_t_derivative_pred=[]
    for x in line_x_axis:
        k_t_derivative_pred.append(__linear__(x,k_linear_args[0],k_linear_args[1]))
    R_2_k_linear=__determination_coefficient__(k_t_derivative,k_t_derivative_pred)
    
   # print(k_linear_args,R_2_k_linear)
   
    k_t_SD,k_linear_SD=__obtainSDfrompcov__(pcov_1),__obtainSDfrompcov__(pcov_2)
    rate_const=k_linear_args[0]*(k_inf-k_0)/c_0
    rate_const_SD=abs(k_linear_SD[0]*(k_inf-k_0)/c_0)
    
    return {
        "file_name":filepath,
        "rate_constant":{"value":rate_const,"standard_deviation":rate_const_SD,"initial_conductivity":k_0,"end_conductivity":k_inf,"initial_concentration":c_0},
            "conductivity_time_calibrating_coefficients":{"values":[pred_args[0],pred_args[1],pred_args[2]],"standard_deviation":k_t_SD,"R_square":R_2_k_t},
            "derivative_conductivity_calibrating_coefficients":{"values":[k_linear_args[0],k_linear_args[1]],"standard_deviation":k_linear_SD,"R_square":R_2_k_linear},
            "calibrated_data":{"time":time,"conductivity":k_t_pred,"derivative":k_t_derivative},
            }
    
def get_activation_energy(filenames,temperatures,k_0,k_inf,c_0,reportname='report'):
    
    rate_constant_data,lnk,T_reciprocal,count=[],[],[],0
    for temp in temperatures:
        T_reciprocal.append(1/temp)
    if len(filenames)==len(temperatures) and len(filenames)==len(k_0) and len(filenames)==len(k_inf) and len(filenames)==len(c_0):
        
        while count<len(filenames):
            data=get_rate_constant(filenames[count],k_0[count],k_inf[count],c_0[count])
            rate_constant_data.append(data)
            lnk.append(math.log(data["rate_constant"]["value"]))
            count+=1
          
        print('Calculating activation energy...')  
        
        calibration_parameters,cov=optimize.curve_fit(__linear__,T_reciprocal,lnk)
        lnk_pred=[]
        for T in T_reciprocal:
            lnk_pred.append(__linear__(T,calibration_parameters[0],calibration_parameters[1]))
        
        R_square_Ea=__determination_coefficient__(lnk,lnk_pred)
        parameters_SD=__obtainSDfrompcov__(cov)
        calibration_parameters=[calibration_parameters[0],calibration_parameters[1]]
        
        activation_energy=-calibration_parameters[0]*molar_gas_constant/1000
        activation_energy_SD=abs(parameters_SD[0]*molar_gas_constant)/1000
        
        print('Activation energy:'+str(activation_energy)+' kJ/mol')
        print('Standard deviation:'+str(activation_energy_SD)+' kJ/mol')
        print('R_square:'+str(R_square_Ea))
        
        plt.scatter(T_reciprocal,lnk,color='black')
        plt.plot(T_reciprocal,lnk_pred,color='blue')
        plt.tick_params(axis='x', labelsize=8)
        plt.tick_params(axis='y', labelsize=8)
        plt.xlabel('1/T')
        plt.ylabel('ln k')
        plt.savefig(reportname+'.png')
        
        reportdata={
            "activation_energy":activation_energy,"standard_deviation":activation_energy_SD,"R_square":R_square_Ea,
                    "lnk_1/T_calibration_coefficients":{"values":calibration_parameters,"standard_deviation":parameters_SD,"1/Temperature":T_reciprocal,"ln_k":lnk},
                    "rate_constant_data":rate_constant_data
                    }
        
        __write_report__(reportname+'.json',reportdata)
        
        print('More data is preserved in '+ reportname +'.json')
        print('Calibration curve is saved in '+ reportname+'.png')
        
        return reportdata
        
    else:
        print("Check your data! Dimensions of arguments in 'get_activation_energy' isn't the same!")
    
    
    
def __obtainSDfrompcov__(pcov):
    
    SD,count=[],0
    while count<len(pcov):
        SD.append(math.sqrt(pcov[count][count]))
        count+=1
    return SD
    
def __determination_coefficient__(y_data,y_pred):
    data_count=len(y_data)
    sum_y,SSres,SStot=0,0,0
    for y in y_data:
        sum_y+=y
    y_mean=sum_y/data_count
    
    count=0
    while count<data_count:
        SSres+=(y_data[count]-y_pred[count])**2
        SStot+=(y_data[count]-y_mean)**2
        count+=1
    
    return (1-SSres/SStot)

def __remove_first_items__(array,count):
    i=0
    while i<count:
        array.pop(i)
        i+=1
        
def __data_pretreatment__(time,k_t):
    count,difference=1,[]
    data_count=len(k_t)
    while count<data_count:
        difference.append(k_t[count]-k_t[count-1])
        count+=1
        
    count_1=0
    while count_1<bad_data_popping_times:
        count,bad_datas=0,bad_data_initial_value
        while count<data_count-1:
            if abs(difference[count])>abs(bad_data_threshold*__medium__(difference)):
                bad_datas.append(count+1)
            count+=1
    
    
        new_time,new_kt,count=[],[],0
        while count<data_count:
            if count in bad_datas:
                pass
            else:
                new_time.append(time[count])
                new_kt.append(k_t[count])
            count+=1
            
        count_1+=1
        
   # print(bad_datas,k_t,difference)
    
    return new_time,new_kt
        
def __medium__(data):
    sorted_data=data.copy()
    sorted_data.sort()
    if len(data)/2-math.floor(len(data)/2)!=0:
        return (sorted_data[int(len(data)/2)]-sorted_data[int(len(data)/2-1)])/2
    else:
        return sorted_data[int((len(data)-1)/2)]

def __write_report__(report_name,dict):
    with open(report_name, 'w', encoding='utf-8') as file:
        json.dump(dict,file,indent=4)
        
get_activation_energy(['25.xlsx','26.xlsx','28.xlsx'],[298.15,299.15,301.15],[2313,2359,2424],[759,780,812],[0.0094,0.0094,0.0094])