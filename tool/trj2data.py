import sys

if len(sys.argv) != 4:
    print("USAGE: python trj2data.py in.lammpstrj OUTPUT_PREFIX GAP")
    sys.exit(-1)

GAP = int(sys.argv[3])

with open(sys.argv[1], 'r') as fp_in:
    counter = 0
    while True:
        line = fp_in.readline()
        if line == "": break
        elif line.startswith("ITEM: NUMBER OF ATOMS"):
            counter += 1
            if (counter % GAP) != 0: continue

            line = fp_in.readline()
            numAtom = int(line.split()[0])
            fp_out = open(f"{sys.argv[2]}{counter}.lammps_data", 'w')
            
            fp_out.write("Comment\n")
            fp_out.write("\n")
            fp_out.write(f"{numAtom} atoms\n")
            fp_out.write(f"3 atom types\n")
            fp_out.write(f"\n")

        elif line.startswith("ITEM: BOX BOUNDS") and ((counter % GAP) == 0):
            line = fp_in.readline()
            splitted = line.split()
            fp_out.write(f"{float(splitted[0])} {float(splitted[1])} xlo xhi\n")

            line = fp_in.readline()
            splitted = line.split()
            fp_out.write(f"{float(splitted[0])} {float(splitted[1])} ylo yhi\n")

            line = fp_in.readline()
            splitted = line.split()
            fp_out.write(f"{float(splitted[0])} {float(splitted[1])} zlo zhi\n")

            fp_out.write(f"\n")

            fp_out.write("Masses\n")
            fp_out.write(f"\n")
            fp_out.write(f"1 6.941\n")
            fp_out.write(f"2 35.453\n")
            fp_out.write(f"3 26.981538\n")
            fp_out.write(f"\n")

        elif line.startswith("ITEM: ATOMS") and ((counter % GAP) == 0):
            fp_out.write(f"Atoms # atomic\n")
            fp_out.write(f"\n")
            for _ in range(numAtom):
                fp_out.write(fp_in.readline())
            fp_out.close()

print("Done! >:D")
