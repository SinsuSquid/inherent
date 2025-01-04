import sys

if len(sys.argv) != 4:
    print("USAGE : python energy.py ENSEMBLE_DIRECTORY LENGTH OUTPUT")
    sys.exit(-1)

ROOT = sys.argv[1]
numTraj = int(sys.argv[2])

fp_out = open(sys.argv[3], 'w')
fp_out.write("#\ttimestep\tPotEng\n")

for i in range(1,numTraj+1):
    with open(f'{ROOT}/DATA/{i}/log.lammps') as fp_in:
        counter = 0
        while True:
            line = fp_in.readline()
            if line == "": break
            elif line.startswith('   Step'):
                if counter == 0:
                    counter += 1
                else:
                    dataline = fp_in.readline()
                    step, time, temp, poteng, kineng, toteng, press, volume, density = dataline.split()

                    toteng = float(toteng)
                    fp_out.write(f"{i}\t{toteng}\n")

fp_out.close()

print("Done >:D!")