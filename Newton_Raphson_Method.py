import numpy as np
import cmath 
j = complex(0,1)

is_slack = False
is_loadbus = False
is_gen = False
initialised = False
continue_as_gen = False

def removing_elements_ifgenbus(J2,J3,J4):
    J2 = 0
    J3 = 0
    J4 = 0
    return J2, J3, J4
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
        B = table[i][4]; #%Succeptance

        #Calculation of Admittance between each line
        Z = R + X
        Admittance = (1/Z) + B

        #Updating the admittance matrix
        Admittance_Array[Bus_from][Bus_to] = Admittance_Array[Bus_from][Bus_to] - Admittance +B
        Admittance_Array[Bus_to][Bus_from] = Admittance_Array[Bus_to][Bus_from] - Admittance +B
        Admittance_Array[Bus_from][Bus_from] = Admittance_Array[Bus_from][Bus_from] + Admittance
        Admittance_Array[Bus_to][Bus_to] = Admittance_Array[Bus_to][Bus_to] + Admittance

   # print(Admittance_Array)
    return Admittance_Array
def initialization_of_voltage(arr_bus):
    for i in range(len(arr_bus)):
        if arr_bus[i][3] == 0:
            arr_bus[i][3] = 1
    return arr_bus
def detect_bus_type(Bus_specification):
    for bus_number in range(len(Bus_specification)):
        if Bus_specification[bus_number][1] == Bus_specification[bus_number][2]:
            Bus_specification[bus_number][5] = "slack"
            print("This is slack\n")
            is_slack = True
        elif Bus_specification[bus_number][3] == Bus_specification[bus_number][4]:
            Bus_specification[bus_number][5] = "pq"
            print("This is laod\n")
            is_loadbus = True

        elif  Bus_specification[bus_number][2] ==  Bus_specification[bus_number][4]:#generator
            Bus_specification[bus_number][5] = "gen"
            print("This is gen\n")
            is_gen = True
def calculate_del_P2(Bus_specification,Y_array, bus_number):
    #p2 calculated for bus_number:
    i = bus_number
    P2_cal = 0

    for j in range(len(Bus_specification)):
        Vi = Bus_specification[i][3]
        Yij = Y_array[i][j]
        Vj = Bus_specification[j][3]
        P2_cal = P2_cal+ abs(Vi)*abs(Vj)*abs(Yij)*np.cos(np.angle(Yij)+np.angle(Vj)-np.angle(Vi))
    P2_specified = Bus_specification[i][1]

    
    del_p = P2_specified - P2_cal
    return del_p
def calculate_del_Q2(Bus_specification, Y_array, bus_number,qmax,qmin):
    #p2 calculated for bus_number:
    i = bus_number
    Q2_cal = 0

    for j in range(len(Bus_specification)):
        Vi = Bus_specification[i][3]
        Yij = Y_array[i][j]
        Vj = Bus_specification[j][3]
        Q2_cal = Q2_cal - abs(Vi)*abs(Vj)*abs(Yij)*np.sin(np.angle(Yij)+np.angle(Vj)-np.angle(Vi))
    cal_val_Q = Q2_cal

    
    Q2_specified = Bus_specification[i][2]


    continue_as_gen = False
    if Bus_specification[bus_number][5] ==  "gen":
        if cal_val_Q > qmin and cal_val_Q < qmax:
            continue_as_gen = True
        elif cal_val_Q < qmin and cal_val_Q < qmax:
            Q2_specified = qmin
        elif cal_val_Q > qmin and cal_val_Q > qmax:
            Q2_specified = qmax
        else:
            print("error in Q")

    del_Q = Q2_specified - Q2_cal

    return del_Q, continue_as_gen
def Jacobian_matrix (jacobian_size,Bus_specification,Y_array,continue_as_gen):
    #here jacobian size = number_of_busses_without_slack
    
    k = 0
    n = 0
    J_matrix = np.zeros((jacobian_size,jacobian_size))
    #print(J_matrix)
    m = int(jacobian_size/2)

    for k in range(len(Bus_specification)):
        for n in range(len(Bus_specification)):
            if n == 0 or k == 0:
                print("jacobian is not for slack bus")
            else:
                if k != n:
                    Vk = abs(Bus_specification[k][3])
                    Vn = abs(Bus_specification[n][3])
                    Ykn = abs(Y_array[k][n])
                    angle_Vk = np.angle(Bus_specification[k][3])
                    angle_Vn = np.angle(Bus_specification[n][3])
                    angle_Ykn = np.angle(Y_array[k][n])

                    J1 = Vk*Vn*Ykn*np.sin(angle_Vk-angle_Vn-angle_Ykn)
                    J2 = Vk*Ykn*np.cos(angle_Vk-angle_Vn-angle_Ykn)
                    J3 = -1*Vk*Vn*Ykn*np.cos(angle_Vk-angle_Vn-angle_Ykn)
                    J4 = Vk*Ykn*np.sin(angle_Vk-angle_Vn-angle_Ykn)
                    if continue_as_gen:
                        J2, J3, J4 = removing_elements_ifgenbus(J2,J3,J4)
                    J_matrix[k-1][n -1] = J1
                    J_matrix[k-1][m +n-1] = J2
                    J_matrix[k+m -1][n-1] = J3
                    J_matrix[k+m -1][m+ n-1] = J4
                    
                    

                elif k==n:
                    J1 = 0

                    for N in range(len(Bus_specification)):
                        if N == k:
                            print("nothing")
                        else:
                            Vk = abs(Bus_specification[k][3])
                            Vn = abs(Bus_specification[N][3])
                            Ykn = abs(Y_array[k][N])
                            angle_Vk = np.angle(Bus_specification[k][3])
                            angle_Vn = np.angle(Bus_specification[N][3])
                            angle_Ykn = np.angle(Y_array[k][N])

                            J1 += Ykn*Vn*np.sin(angle_Vk-angle_Vn-angle_Ykn)
                    
                    J1 = J1 * -Vk
                    
                    J2 = 0

                    for N in range(len(Bus_specification)):

                        Vk = abs(Bus_specification[k][3])
                        Vn = abs(Bus_specification[N][3])
                        Ykn = abs(Y_array[k][N])
                        angle_Vk = np.angle(Bus_specification[k][3])
                        angle_Vn = np.angle(Bus_specification[N][3])
                        angle_Ykn = np.angle(Y_array[k][N])

                        J2 += Ykn*Vn*np.cos(angle_Vk-angle_Vn-angle_Ykn)
                    
                    J2 = J2 + Vk*abs(Y_array[k][k])*np.cos(np.angle(Y_array[k][k]))

                    J4 = 0

                    for N in range(len(Bus_specification)):

                        Vk = abs(Bus_specification[k][3])
                        Vn = abs(Bus_specification[N][3])
                        Ykn = abs(Y_array[k][N])
                        angle_Vk = np.angle(Bus_specification[k][3])
                        angle_Vn = np.angle(Bus_specification[N][3])
                        angle_Ykn = np.angle(Y_array[k][N])

                        J4 += Ykn*Vn*np.sin(angle_Vk-angle_Vn-angle_Ykn)
                    
                    J4 = J4 - Vk*abs(Y_array[k][k])*np.sin(np.angle(Y_array[k][k]))

                    J3 = 0

                    for N in range(len(Bus_specification)):
                        if N == k:
                            print("nothing")
                        else:
                            Vk = abs(Bus_specification[k][3])
                            Vn = abs(Bus_specification[N][3])
                            Ykn = abs(Y_array[k][N])
                            angle_Vk = np.angle(Bus_specification[k][3])
                            angle_Vn = np.angle(Bus_specification[N][3])
                            angle_Ykn = np.angle(Y_array[k][N])

                            J3 += Ykn*Vn*np.cos(angle_Vk-angle_Vn-angle_Ykn)
                    
                    J3 = J3 * Vk
                    if continue_as_gen:
                        J2, J3, J4 = removing_elements_ifgenbus(J2,J3,J4)

                    J_matrix[k-1][n -1] = J1
                    J_matrix[k-1][m +n-1] = J2
                    J_matrix[k+m -1][n-1] = J3
                    J_matrix[k+m -1][m+ n-1] = J4
                
    return J_matrix





#Input From User without 'j' % Total succeptance per line is taken as input (For partb,
#succeptance has to be doubled before inputing into program)
table = [

 [1, 2 ,2.8, -9.6, 2.799-9.6001*j] #[busfrom, bus to, resistance, reactance, succeptance (total)]
 #[1 ,3, 0.02, 0.06, 0],
 #[2 ,3 ,0.06 ,0.18 ,0]#,
 #[2 ,4 ,0.08 ,0.24 ,0.05],
 #[2 ,5 ,0.02 ,0.06 ,0.02],
 #[3 ,4 ,0.01 ,0.04 ,0.01],
 #[4 ,5 ,0.03 ,0.10 ,0.04]
]

Y_array = find_admittance_array (table)


# for below, enter 0 for the values that are not given
Bus_specification = [

    [1, 0, 0, 1.0, 0, "slack"],#[bus code, P , Q, V, del(in degrees)]
    [2, -1.5, -0.5, 0, 0, "pq"],
    #[3, 0.6, 0.25, 0, 0, "pq"]#,
    #[1, 0, 0, 0, 0, "pq"],
    #[1, 0, 0, 0, 0, "pq"],
    #[1, 0, 0, 0, 0, "pq"],
    #[1, 0, 0, 0, 0, "pq"],
    #[1, 0, 0, 0, 0, "pq"]
]
number_of_busses_without_slack = len(Bus_specification)

qmin = 0
qmax = 0.35
num_iterations = 1

jacobian_size = 0

jacobian_size = detect_bus_type(Bus_specification)


del_p_q_array = np.zeros(2*(number_of_busses_without_slack-1))

for n in range(num_iterations):
    for bus_number in range(len(Bus_specification)):#Slack
        if Bus_specification[bus_number][5] == "slack" :
            slack_bus_num = bus_number
        elif Bus_specification[bus_number][5] == "pq":#PQ
            load_bus = bus_number
            print(load_bus)
            if not(initialised):
                Bus_specification = initialization_of_voltage(Bus_specification) #Initializing voltage
                
            else:
                Bus_specification = Bus_specification
                
            print("hellp1" +str(number_of_busses_without_slack + load_bus - 2))
            del_p_q_array[load_bus - 1] = calculate_del_P2(Bus_specification,Y_array, load_bus)
            del_p_q_array[number_of_busses_without_slack + load_bus - 2],continue_dummy_variable= calculate_del_Q2(Bus_specification, Y_array, load_bus,qmax,qmin)
        
            
        elif  Bus_specification[bus_number][5] ==  "gen":#generator
            gen_bus = bus_number
            print("hellp" +str(number_of_busses_without_slack + gen_bus - 2))
            previous_bus_voltage = Bus_specification[gen_bus][3]
            q_temp, continue_as_gen = calculate_del_Q2(Bus_specification, Y_array, gen_bus,qmax,qmin)
            del_p_q_array[gen_bus - 1] = calculate_del_P2(Bus_specification,Y_array, gen_bus)
            del_p_q_array[number_of_busses_without_slack + gen_bus - 2]= q_temp

            if continue_as_gen:
                print(continue_as_gen)
                Bus_specification[gen_bus][3] = previous_bus_voltage
            else:
                Bus_specification[gen_bus][3] = 1

    Final_J_Matrix = Jacobian_matrix (number_of_busses_without_slack,Bus_specification,Y_array,continue_as_gen)
    if Bus_specification[bus_number][5] ==  "gen":
        Bus_specification[gen_bus][3] = previous_bus_voltage

print(Final_J_Matrix)
print(del_p_q_array)
if not(continue_as_gen):
    del_V_del =np.matmul(np.linalg.inv(Final_J_Matrix),del_p_q_array)
else:
    del_del = del_p_q_array[0]/Final_J_Matrix[0][0]

print(del_V_del)


