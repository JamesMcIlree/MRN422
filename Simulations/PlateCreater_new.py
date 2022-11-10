import math
import random
import statistics
import gmsh
import os
import matplotlib.pyplot as plt
import numpy as np

def get_constants(): # Enter Constants
    number_of_plates = 1

    plate_width = 200/1000 # m
    plate_length = 200/1000 # m
    plate_thickness = 1.95/1000 # m

    mesh_plate = 5e-3 # Fineness of mesh around plate
    mesh_plate_crack = 5e-3 # Fineness of mesh around plate with crack
    mesh_crack = 1e-3 # Fineness of mesh around crack

    upper_crack_limit = 60/1000 # m
    lower_crack_limit = 40/1000 # m

    return number_of_plates, plate_width, plate_length, plate_thickness, mesh_plate, mesh_plate_crack, mesh_crack, lower_crack_limit, upper_crack_limit

def position_on_plate(no): # Determines the quater the crack will generate on
    pos = 0
    count = 0
    positions = []

    while count < no:
        positions.append(pos)
        pos = pos + 1
        count = count +1

        if pos > 4:
            pos = 0    

    return positions

def randomize_crack_variables(lcl, ucl):
    crack_length = random.uniform(lcl, ucl)
    crack_width = random.uniform(lcl*0.1,ucl*0.2)
    crack_angle = random.uniform(0, 2*np.pi)
    return round(crack_length, 10), round(crack_width, 10), round(crack_angle, 10)

def random_crack_position(position, cl, pw, pl): # Creates a random starting location for crack
    if position == 1: # Top right
        box_x1 = pw/2 - 2*cl
        box_x2 = 0 + 2*cl

        box_y1 = pl/2 - 2*cl
        box_y2 = 0 + 2*cl

    elif position == 2: # Top left
        box_x1 = -pw/2 + 2*cl
        box_x2 = 0 - 2*cl

        box_y1 = pl/2 - 2*cl
        box_y2 = 0 + 2*cl

    elif position == 3: # bottom left
        box_x1 = -pw/2 + 2*cl
        box_x2 = 0 - 2*cl

        box_y1 = -pl/2 + 2*cl
        box_y2 = 0 - 2*cl

    elif position == 4: # bottom right
        box_x1 = pw/2 - 2*cl
        box_x2 = 0 + 2*cl

        box_y1 = -pl/2 + 2*cl
        box_y2 = 0 - 2*cl

    else: 
        box_x1 = 0
        box_x2 = 0

        box_y1 = 0
        box_y2 = 0

    crack_x = random.uniform(box_x1, box_x2)
    crack_y = random.uniform(box_y2, box_y1)
    
    return round(crack_x,3), round(crack_y,3)

def get_crack_coords(cl, cw, ctheta, pos, pw, pl): # Creates coordinates of crack
    check = True
    while check:
        crack_coord_x1, crack_coord_y1 = random_crack_position(pos, cl, pw, pl)

        crack_coord_x2 = crack_coord_x1 + cw*math.cos(ctheta)
        crack_coord_x3 = crack_coord_x1 - cw*math.cos(ctheta)

        crack_coord_y2 = crack_coord_y1 - cw*math.sin(ctheta)
        crack_coord_y3 = crack_coord_y1 + cw*math.sin(ctheta)

        crack_coord_x4 = crack_coord_x1 + cl*math.cos(ctheta)
        crack_coord_y4 = crack_coord_y1 + cl*math.sin(ctheta)

        crack_coord_x5 = crack_coord_x4 + cw*math.cos(ctheta)
        crack_coord_x6 = crack_coord_x4 - cw*math.cos(ctheta)

        crack_coord_y5 = crack_coord_y4 - cw*math.sin(ctheta)
        crack_coord_y6 = crack_coord_y4 + cw*math.sin(ctheta)

        if np.abs(crack_coord_x1) < pw/2 and np.abs(crack_coord_x2) < pw/2 and np.abs(crack_coord_x3) < pw/2 and np.abs(crack_coord_x4) < pw/2 and np.abs(crack_coord_x5) < pw/2 and np.abs(crack_coord_x6) < pw/2:
            if np.abs(crack_coord_y1) < pl/2 and np.abs(crack_coord_y2) < pl/2 and np.abs(crack_coord_y3) < pl/2 and np.abs(crack_coord_y4) < pl/2 and np.abs(crack_coord_y5) < pl/2 and np.abs(crack_coord_y6) < pl/2:
                check = False

    return crack_coord_x1, crack_coord_x2, crack_coord_x3, crack_coord_x4, crack_coord_x5, crack_coord_x6, crack_coord_y1, crack_coord_y2, crack_coord_y3, crack_coord_y4, crack_coord_y5, crack_coord_y6
    
def make_INP_file(no, pw, pl, pt, mp, mpc, mc, lcl, ucl): # makes a INP from gmsh for ccx
    positions = position_on_plate(no)
    crack_det = []

    counter = 1
    while counter <= no:
        cl, cw, ctheta = randomize_crack_variables(lcl,ucl)
        ccx1, ccx2, ccx3, ccx4, ccx5, ccx6, ccy1, ccy2, ccy3, ccy4, ccy5, ccy6 = get_crack_coords(cl, cw, ctheta, positions[counter - 1], pw, pl)
        name = "plate" + str(counter)
        plate_details(cl, cw, ctheta, counter)

        if positions[counter - 1] != 0:
            gmsh.initialize()
            gmsh.model.add(name)

            gmsh.model.geo.addPoint(pl/2,pw/2,0,mpc,1) # Plate corners
            gmsh.model.geo.addPoint(-pl/2,pw/2,0,mpc,2)
            gmsh.model.geo.addPoint(-pl/2,-pw/2,0,mpc,3)
            gmsh.model.geo.addPoint(pl/2,-pw/2,0,mpc,4)

            gmsh.model.geo.addPoint(ccx1,ccy1,0,mc,5) # Crack points
            gmsh.model.geo.addPoint(ccx2,ccy2,0,mc,6)
            gmsh.model.geo.addPoint(ccx3,ccy3,0,mc,7)

            gmsh.model.geo.addPoint(ccx4,ccy4,0,mc,8)
            gmsh.model.geo.addPoint(ccx5,ccy5,0,mc,9)
            gmsh.model.geo.addPoint(ccx6,ccy6,0,mc,10)

            gmsh.model.geo.addLine(1,2,11) # Plate edges
            gmsh.model.geo.addLine(2,3,12)
            gmsh.model.geo.addLine(3,4,13)
            gmsh.model.geo.addLine(4,1,14)

            gmsh.model.geo.addLine(7,6,15) # Crack lines
            gmsh.model.geo.addLine(10,9,16)

            gmsh.model.geo.addLine(7,10,17)
            gmsh.model.geo.addLine(6,9,18)

            gmsh.model.geo.addCurveLoop([11,12,13,14], 19)
            gmsh.model.geo.addCurveLoop([15,-17,-16,18], 20)

            s = gmsh.model.geo.addPlaneSurface([19,20], 21)
            e = gmsh.model.geo.extrude([(2,s)], 0, 0, pt)
            p = gmsh.model.geo.addPhysicalGroup(3, [1], 1)

            gmsh.model.geo.addPoint(0,0,0,mpc,10000) # Shaker point

            gmsh.model.geo.addPoint(-pw/4,pl/4,pt,mpc,10001) # Accelerometer Points
            gmsh.model.geo.addPoint(pw/4,pl/4,pt,mpc,10002)
            gmsh.model.geo.addPoint(0,0,pt,mpc,10003)
            gmsh.model.geo.addPoint(-pw/4,-pl/4,pt,mpc,10004)
            gmsh.model.geo.addPoint(pw/4,-pl/4,pt,mpc,10005)  

            '''gmsh.model.geo.synchronize()

            gmsh.model.geo.mesh.setTransfiniteSurface(21)
            gmsh.model.geo.mesh.setTransfiniteSurface(30)
            gmsh.model.geo.mesh.setTransfiniteSurface(34)
            gmsh.model.geo.mesh.setTransfiniteSurface(38)
            gmsh.model.geo.mesh.setTransfiniteSurface(42)
            gmsh.model.geo.mesh.setTransfiniteSurface(43)

            gmsh.model.geo.synchronize()

            gmsh.model.geo.mesh.setTransfiniteVolume(1)

            gmsh.model.geo.mesh.setRecombine(2, 21)
            gmsh.model.geo.mesh.setRecombine(2, 34)
            gmsh.model.geo.mesh.setRecombine(2, 38)
            gmsh.model.geo.mesh.setRecombine(2, 42)
            gmsh.model.geo.mesh.setRecombine(2, 46)
            gmsh.model.geo.mesh.setRecombine(2, 63)
            gmsh.model.geo.mesh.setRecombine(3, 1)'''

            gmsh.model.geo.synchronize()

            gmsh.model.mesh.embed(0,[10000],2,21)

            gmsh.model.mesh.embed(0,[10001],2,63)
            gmsh.model.mesh.embed(0,[10002],2,63)
            gmsh.model.mesh.embed(0,[10003],2,63)
            gmsh.model.mesh.embed(0,[10004],2,63)
            gmsh.model.mesh.embed(0,[10005],2,63)

            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.setOrder(2)
            gmsh.write("plate" + str(counter) + ".inp")
            #gmsh.fltk.run()
            gmsh.finalize()
            crack_det.append(1)

        else:
            gmsh.initialize()
            gmsh.model.add(name)

            gmsh.model.geo.addPoint(pl/2,pw/2,0,mp,1) # Plate corners
            gmsh.model.geo.addPoint(-pl/2,pw/2,0,mp,2)
            gmsh.model.geo.addPoint(-pl/2,-pw/2,0,mp,3)
            gmsh.model.geo.addPoint(pl/2,-pw/2,0,mp,4)

            gmsh.model.geo.addLine(1,2,11) # Plate edges
            gmsh.model.geo.addLine(2,3,12)
            gmsh.model.geo.addLine(3,4,13)
            gmsh.model.geo.addLine(4,1,14)

            gmsh.model.geo.addCurveLoop([11,12,13,14], 19)

            s = gmsh.model.geo.addPlaneSurface([19], 21)
            e = gmsh.model.geo.extrude([(2,s)], 0, 0, pt)
            p = gmsh.model.geo.addPhysicalGroup(3, [1], 1)

            gmsh.model.geo.addPoint(0,0,0,mp,10000) # Shaker point

            gmsh.model.geo.addPoint(-pw/4,pl/4,pt,mp,10001) # Accelerometer Points
            gmsh.model.geo.addPoint(pw/4,pl/4,pt,mp,10002)
            gmsh.model.geo.addPoint(0,0,pt,mp,10003)
            gmsh.model.geo.addPoint(-pw/4,-pl/4,pt,mp,10004)
            gmsh.model.geo.addPoint(pw/4,-pl/4,pt,mp,10005)

            '''gmsh.model.geo.synchronize()

            gmsh.model.geo.mesh.setTransfiniteSurface(21, "Left")
            gmsh.model.geo.mesh.setTransfiniteSurface(30, "Left")
            gmsh.model.geo.mesh.setTransfiniteSurface(34, "Left")
            gmsh.model.geo.mesh.setTransfiniteSurface(38, "Left")
            gmsh.model.geo.mesh.setTransfiniteSurface(42, "Left")
            gmsh.model.geo.mesh.setTransfiniteSurface(43, "Left")

            gmsh.model.geo.synchronize()

            gmsh.model.geo.mesh.setTransfiniteVolume(1)

            gmsh.model.geo.mesh.setRecombine(2, 21)
            gmsh.model.geo.mesh.setRecombine(2, 30)
            gmsh.model.geo.mesh.setRecombine(2, 34)
            gmsh.model.geo.mesh.setRecombine(2, 38)
            gmsh.model.geo.mesh.setRecombine(2, 42)
            gmsh.model.geo.mesh.setRecombine(2, 43)
            gmsh.model.geo.mesh.setRecombine(3, 1)'''

            gmsh.model.geo.synchronize()

            gmsh.model.mesh.embed(0,[10000],2,21)

            gmsh.model.mesh.embed(0,[10001],2,43)
            gmsh.model.mesh.embed(0,[10002],2,43)
            gmsh.model.mesh.embed(0,[10003],2,43)
            gmsh.model.mesh.embed(0,[10004],2,43)
            gmsh.model.mesh.embed(0,[10005],2,43)

            gmsh.model.geo.synchronize()

            gmsh.model.mesh.generate(3)
            gmsh.model.mesh.setOrder(2)
            gmsh.write("plate" + str(counter) + ".inp")
            #gmsh.fltk.run()
            gmsh.finalize()
            crack_det.append(0)

        counter = counter + 1
    
    return crack_det

def edit_INP(no, crack_det): # Edits INP for ccx
    counter = 1

    while counter <= no:
        file_name = "plate" + str(counter) + ".inp"
        plate = open(file_name, "a")

        line1 = "*MATERIAL, Name=steel\n*ELASTIC\n"
        line2 = "209e9, 0.292\n"
        line3 = "*DENSITY\n"
        line4 = "7800\n"
        line5 = "*SOLID SECTION,ELSET=Volume1,MATERIAL=steel\n"
        line6 = "*STEP\n"
        line7 = "*FREQUENCY\n"
        line8 = "15,1,1000\n"
        line9 = "*NODE FILE\n"
        line10 = "U\n"
        line11 = "*EL FILE\n"
        line12 = "E, S\n"
        line13 = "*END STEP"

        lines = [line1,line2,line3,line4,line5,line6,line7,line8,line9,line10,line11,line12,line13]

        for line in lines:
            plate.write(line)

        plate.close()
        counter = counter + 1

    return no

def plate_details(cl, cw, ctheta, plate):
    name = "details.txt"
    file = open(name,"a")
    file.write("Plate " + str(plate) + "\n")
    file.write("Crack Length = " + str(cl*1000) + "mm\n")
    file.write("Crack Width = " + str(cw*1000) + "mm\n")
    file.write("Crack Angle = " + str(ctheta*180/np.pi) + "\n\n")
    file.close()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def solver(no, crack_det):
    counter = 1

    while counter <= no:
        command1 = "ccx plate" + str(counter)
        command2 = "del /f plate" + str(counter) + ".CVG"
        command3 = "del /f plate" + str(counter) + ".DAT"
        command4 = "del /f plate" + str(counter) + ".EIG"
        command5 = "del /f plate" + str(counter) + ".STA"
        command6 = "del /f spooles.OUT"

        os.system(command1)
        os.system(command2)
        os.system(command3)
        os.system(command4)
        os.system(command5)
        os.system(command6)
        get_node_data(counter, crack_det[counter-1])
        counter = counter + 1

def get_node_data(no, crack_det):
    if crack_det == 0:
        make_command_file_no_crack(no)
        command1 = "cgx -bg command" + str(no) + ".fbd"
        command2 = "del /f graph_0.gnu"
        command3 = "del /f graph_0.OUT"
        command4 = "del /f graph_0.OUT2"
        command5 = "del /f graph_0.PNG"
        command6 = "del /f command" + str(no) + ".fbd"
        command7 = "del /f plate" + str(no) + ".frd"
        command8 = "del /f plate" + str(no) + ".inp"

        os.system(command1)
        os.system(command2)
        os.system(command3)
        os.system(command4)
        os.system(command5)
        os.system(command6)
        os.system(command7)
        os.system(command8)

    else:
        make_command_file_crack(no)
        command1 = "cgx -bg command" + str(no) + ".fbd"
        command2 = "del /f graph_0.gnu"
        command3 = "del /f graph_0.OUT"
        command4 = "del /f graph_0.OUT2"
        command5 = "del /f graph_0.PNG"
        command6 = "del /f command" + str(no) + ".fbd"
        command7 = "del /f plate" + str(no) + ".frd"
        command8 = "del /f plate" + str(no) + ".inp"

        os.system(command1)
        os.system(command2)
        os.system(command3)
        os.system(command4)
        os.system(command5)
        os.system(command6)
        os.system(command7)
        os.system(command8)

def make_command_file_crack(number):
    file = open("command" + str(number) + ".fbd","w")

    line1 = "read plate" + str(number) + ".frd\n"
    line2 = "seta set" + str(number) + " n 18 19 20 21 22\n"
    line3 = "graph set" + str(number) + " time DISP ALL"

    lines = [line1, line2, line3]

    for line in lines:
        file.write(line)

    file.close()

def make_command_file_no_crack(number):
    file = open("command" + str(number) + ".fbd","w")

    line1 = "read plate" + str(number) + ".frd\n"
    line2 = "seta set" + str(number) + " n 10 11 12 13 14\n"
    line3 = "graph set" + str(number) + " time DISP ALL"

    lines = [line1, line2, line3]

    for line in lines:
        file.write(line)

    file.close()

#---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

def postprocessing(no):
    counter = 1
    plates = []
    while counter <= no:
        counter2 = 1
        plate_pos = []
        plate_vel = []
        plate_acc = []
        time, results = read_results(counter)
        for node_pos in results:
            vel, dx = make_velocity_function(node_pos, time)
            acc, dx2 = make_acceleration_function(vel, dx)
            display_results(time, node_pos, counter, counter2)
            plate_pos.append(node_pos)
            plate_vel.append(vel)
            plate_acc.append(acc)
            counter2 = counter2 + 1
        counter = counter + 1
        plates.append([plate_pos, plate_vel, plate_acc])
    make_stats(plates)
        
def read_results(no):
    lines = []
    file = open("graph_set" + str(no) + "_DISP_ALL.OUT","r")
    for line in file:
        lines.append(line)
    file.close()
    time, results = seperate_results(lines)
    return time, results

def seperate_results(lines):
    time = []
    node1 = []
    node2 = []
    node3 = []
    node4 = []
    node5 = []
    node6 = []
    node7 = []
    node8 = []
    node9 = []
    counter = 0

    for line in lines:
        sep_line = line.split(";")
        if len(sep_line) > 8:
            if counter != 0:
                time.append(float(sep_line[0].strip()))
                node1.append(float(sep_line[1].strip()))
                node2.append(float(sep_line[2].strip()))
                node3.append(float(sep_line[3].strip()))
                node4.append(float(sep_line[4].strip()))
                node5.append(float(sep_line[5].strip()))
                node6.append(float(sep_line[6].strip()))
                node7.append(float(sep_line[7].strip()))
                node8.append(float(sep_line[8].strip()))
                node9.append(float(sep_line[9].strip()))
            counter = 1

    nodes = [node1,node2,node3,node4,node5,node6,node7,node8,node9]
    return time, nodes

def display_results(time, result, plate, node):
    y = result
    x = time
    plt.title("Plate " + str(plate) + " Node " + str(node), fontsize=18)
    plt.xlabel("Time [s]", fontsize=18)
    plt.ylabel("Acceleration [m/s^2]", fontsize=18)
    plt.plot(x,y,"b-")
    plt.show()

def make_stats(plates):
    name = "stats.txt"
    file = open(name,"w")
    count1 = 1
    
    for plate in plates:
        file.write("Plate" + str(count1) + "\n")
        count2 = 1

        for var in plate:
            count3 = 1
            if count2 == 1:
                file.write("Position\n")
            elif count2 == 2:
                file.write("Velocity\n")
            else: 
                file.write("Acceleration\n")
            count2 = count2 + 1

            file.write("Node; sum; mean; median; mode; max, min; range\n")

            for node in var:
                Sum = str(sum(node)) + "; "
                mean = str(statistics.mean(node)) + "; "
                median = str(statistics.median(node)) + ";"
                mode = str(statistics.mode(node)) + "; "
                Max = str(max(node)) + "; "
                Min = str(min(node)) + "; "
                range = str(max(node) - min(node)) + "\n"
                file.write(str(count3) + "; " + Sum + mean + median + mode + Max + Min + range)
                count3 = count3 + 1
            
            file.write("\n")
        file.write("\n")
        count1 = count1 + 1

    file.close()

def make_velocity_function(FRF, time):
    dx = np.diff(time)
    velocity_function = np.diff(FRF)/dx
    time = time[:-1]
    return velocity_function, time

def make_acceleration_function(vel, time):
    dx = np.diff(time)
    acceleration_function = np.diff(vel)/dx
    time = time[:-1]
    return acceleration_function,time

def main():
    no, pw, pl, pt, mp, mpc, mc, lcl, ucl = get_constants()
    crack_det = make_INP_file(no, pw, pl, pt, mp, mpc, mc, lcl, ucl)
    edit_INP(no, crack_det)
    #solver(no, crack_det)
    #postprocessing(no)
    
main()