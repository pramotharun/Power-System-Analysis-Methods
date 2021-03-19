import numpy as np 
import cmath 
j = complex(0,1)
def find_max(arr_table):
    old_max = 0
    for i in range(len(arr_table)):
        a = arr_table[i][0]
        b = arr_table[i][1]

        new_max = max(a,b,old_max)
        old_max = new_max
    return old_max
def find_admittance_array (table):
    #To find the total number of busses from the 2 columns
    max_value = find_max(table) # %Gives total number of busses

    #print(max_value)
    #Initialising admittance array with the number of busses
    Admittance_Array = np.zeros((max_value,max_value),dtype = "complex_")

    j = complex(0,1)

    #For loop interating through the different admittances in each row
    for i in range(len(table)):
        # Taking from input:
        Bus_from = table[i][0] - 1
        Bus_to = table[i][1] - 1
        R = table[i][2] #Resistance
        X = table[i][3]* j;#% Reactance
        B = table[i][4]* j/2; #%Succeptance

        #Calculation of Admittance between each line
        Z = R + X
        Admittance = (1/Z) + B

        #Updating the admittance matrix
        Admittance_Array[Bus_from][Bus_to] = Admittance_Array[Bus_from][Bus_to] - Admittance +B
        Admittance_Array[Bus_to][Bus_from] = Admittance_Array[Bus_to][Bus_from] - Admittance +B
        Admittance_Array[Bus_from][Bus_from] = Admittance_Array[Bus_from][Bus_from] + Admittance
        Admittance_Array[Bus_to][Bus_to] = Admittance_Array[Bus_to][Bus_to] + Admittance

    return Admittance_Array
def angle_radian(arr_bus):
    for i in range(len(arr_bus)):
        arr_bus[i][4] = arr_bus[i][4] *3.14/180
    return arr_bus
def initialization_of_voltage(arr_bus):
    for i in range(len(arr_bus)):
        if arr_bus[i][3] == 0:
            arr_bus[i][3] = 1
    return arr_bus
def voltage_admittance_multiplication(bus_number,Bus_specification,final_admittance_matrix):
    sum = 0
    for j in range(len(final_admittance_matrix)):
        if j == bus_number:
            j+= 1
        else:
            sum += final_admittance_matrix[bus_number][j]*Bus_specification[j][3]
    return sum
def calculate_Q(qmin,qmax,bus_number,Bus_specification,final_admittance_matrix, table):
    #cal_val = np.conj(Bus_specification[bus_number][3])* voltage_admittance_multiplication(bus_number,Bus_specification,final_admittance_matrix,True)
    #cal_val_Q = -1*cal_val.imag
    #print("this is calculated Q               "+str(cal_val_Q))
    cal_val_Q = 0
    for i in range(find_max(table)):
        Fp = np.imag(Bus_specification[bus_number][3])
        Ep = np.real(Bus_specification[bus_number][3])
        Fq = np.imag(Bus_specification[i][3])
        Eq = np.real(Bus_specification[i][3])
        Gpq = np.real(final_admittance_matrix[bus_number][i])
        Bpq = np.imag(final_admittance_matrix[bus_number][i])
        cal_val_Q = cal_val_Q + Fp*(Eq*Gpq + Fp*Bpq) - Ep*(Fp*Gpq + Eq*Bpq)
    
    #print("this is calculated Q               "+str(cal_val_Q))
    if cal_val_Q > qmin and cal_val_Q < qmax:
        return cal_val_Q, cal_val_Q
    elif cal_val_Q < qmin and cal_val_Q < qmax:
        return qmin, cal_val_Q
    elif cal_val_Q > qmin and cal_val_Q > qmax:
        return qmax, cal_val_Q
    else:
        print("error in Q")




#Input From User without 'j' % Total succeptance per line is taken as input (For partb,
#succeptance has to be doubled before inputing into program)
table = [

 [1, 2 ,0, 0.020, 0], #[busfrom, bus to, resistance, reactance, succeptance (total)]
 [1 ,3, 0, 0.045, 0],
 [2 ,3 ,0 ,0.020 ,0]
 #[2 ,4 ,0.08 ,0.24 ,0.05],
 #[2 ,5 ,0.02 ,0.06 ,0.02],
 #[3 ,4 ,0.01 ,0.04 ,0.01],
 #[4 ,5 ,0.03 ,0.10 ,0.04]
]

# for below, enter 0 for the values that are not given
Bus_specification = [

    [1, 0, 0, 1.025, 0, "slack"],#[bus code, P , Q, V, del(in degrees)]
    [2, -0.78, -0.38, 0, 0, "pq"],
    [3, 0.58, 0, 1.03, 0, "pq"]
    #[1, 0, 0, 0, 0, "pq"],
    #[1, 0, 0, 0, 0, "pq"],
    #[1, 0, 0, 0, 0, "pq"],
    #[1, 0, 0, 0, 0, "pq"],
    #[1, 0, 0, 0, 0, "pq"]
]

alpha = 1.1 # convergance factor
num_iterations = 1 # Number of iterations
qmin = -20000
qmax = 100000

Bus_specification = angle_radian(Bus_specification) # If we have to convert to radians
initialised  = False
final_admittance_matrix = find_admittance_array(table)
is_gen = False
slack_bus = False
is_loadbus = False

for bus_number in range(len(Bus_specification)):#Slack
    if Bus_specification[bus_number][1] == Bus_specification[bus_number][2]:
        Bus_specification[bus_number][5] = "slack"
        bus_number+=1
        print("entered ere")
    elif Bus_specification[bus_number][3] == Bus_specification[bus_number][4]:
        Bus_specification[bus_number][5] = "pq"
        bus_number+=1
        print("entered ere1")
    elif  Bus_specification[bus_number][2] ==  Bus_specification[bus_number][4]:#generator
        Bus_specification[bus_number][5] = "gen"
        bus_number+=1
        print("entered ere2")
    
    




for n in range(num_iterations):
    for bus_number in range(len(Bus_specification)):#Slack
        if Bus_specification[bus_number][1] == Bus_specification[bus_number][2]:
            print("entered slack")
            #print("Slack bus: "+str(bus_number))
            slack_bus  = True
            slack_bus_num = bus_number
           # bus_number+=1
        elif Bus_specification[bus_number][3] == Bus_specification[bus_number][4]:#PQ
            print("entered Load PQ")
            is_loadbus = True
            load_bus = bus_number
            previous_voltage = Bus_specification[bus_number][3]


            Bus_specification = initialization_of_voltage(Bus_specification) #Initializing voltage

            
            Bus_specification[bus_number][3] = (1*( ((Bus_specification[bus_number][1] - j*Bus_specification[bus_number][2])/np.conj(Bus_specification[bus_number][3]) ) - voltage_admittance_multiplication(bus_number,Bus_specification,final_admittance_matrix)))/final_admittance_matrix[bus_number][bus_number]
            Bus_specification[bus_number][3] = previous_voltage + alpha*(Bus_specification[bus_number][3]- previous_voltage)
            #bus_number+=1
        elif  Bus_specification[bus_number][2] ==  Bus_specification[bus_number][4]:#generator
            gen_bus = bus_number
            is_gen = True
            Q_value, actual = calculate_Q(qmin,qmax,bus_number,Bus_specification,final_admittance_matrix,table)
            Bus_specification[bus_number][2] = actual
            #print(Bus_specification)
    
            if Q_value==qmin or Q_value==qmax:
                previous_voltage = Bus_specification[bus_number][3]
                gen_voltage = ((((Bus_specification[bus_number][1] - j*Q_value)/np.conj(1+j*0) ) - voltage_admittance_multiplication(bus_number,Bus_specification,final_admittance_matrix)))/final_admittance_matrix[bus_number][bus_number]

                Bus_specification[bus_number][3] = 1*(np.cos(np.angle(gen_voltage))+j*np.sin(np.angle(gen_voltage)))
                Bus_specification[bus_number][3] = previous_voltage + alpha*(Bus_specification[bus_number][3]- previous_voltage)
                
            else:
                previous_voltage = Bus_specification[bus_number][3]

                Bus_specification[bus_number][3] = ((((Bus_specification[bus_number][1] - j*Bus_specification[bus_number][2])/np.conj(Bus_specification[bus_number][3]) ) - voltage_admittance_multiplication(bus_number,Bus_specification,final_admittance_matrix) ))/final_admittance_matrix[bus_number][bus_number]
                angle_new = np.angle(Bus_specification[bus_number][3])
                Bus_specification[bus_number][3] = abs(previous_voltage)*(np.cos(angle_new)+j*np.sin(angle_new))
                
                Bus_specification[bus_number][3] = previous_voltage + alpha*(Bus_specification[bus_number][3]- previous_voltage)

        

                

    if is_gen:
        gen_bus = gen_bus
    elif is_loadbus:
        gen_bus = load_bus

    if is_gen or is_loadbus:

        Vi = Bus_specification[slack_bus_num][3] #slack bus voltage
        
        Vj = Bus_specification[gen_bus][3] #generator bus voltage
        # if there is load bus change to: load_bus
        Yij = -1*final_admittance_matrix[slack_bus_num][gen_bus]
        Yji = -1*final_admittance_matrix[gen_bus][slack_bus_num]

        #### shouldnt this Ysi succeptance be the total succeptance that is connected to ground in that bus?

        Ysi = 0
        Ysj = 0

        for i in range(len(table)):
            if table[i][0] == slack_bus_num+1 or table[i][1] == slack_bus_num+1:
                Ysi += table[i][4]/2

            if table[i][0] == (gen_bus+1) or table[i][1] == (gen_bus+1):
                Ysj += table[i][4]/2


        

        Sij = Vi*((np.conj(Vi)-np.conj(Vj))*np.conj(Yij)) + (abs(Vi)**2)*np.conj(Ysi*j) #here p is the bus from
        Sji = Vj*((np.conj(Vj)-np.conj(Vi))*np.conj(Yji)) + (abs(Vj)**2)*np.conj(Ysj*j) #here p is the bus from
        print("s12" +str(Sij))
        print("s12" +str(Sji))
        P_ADD_Q = Sij + Sji

        Bus_specification[slack_bus_num][1] = np.real(P_ADD_Q)
        Bus_specification[slack_bus_num][2] = np.imag(P_ADD_Q)
    
    


print("____________________________________________________________\n")

print("Final bus specifications NOTE UNITs: PU\n")
output_array = ["Bus Number", "P", "Q", "Voltage", "delta"]
print(str(output_array) + "\n")
for i in range(len(Bus_specification)):
    print(str(Bus_specification[i]) + "\n")

print("____________________________________________________________\n")

print(final_admittance_matrix)

