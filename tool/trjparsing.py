import sys

if len(sys.argv) != 3:
    print("USAGE: python trjparsing.py in.lammpstrj out.lammpstrj")
    sys.exit(-1)

fp_out = open(sys.argv[2], 'w')

with open(sys.argv[1], 'r') as fp_in:
    counter = -1
    while True:
        line = fp_in.readline()
        if line == "": break
        elif line.startswith("ITEM: TIMESTEP"):
            counter += 1
            if (counter % 10) != 0:
                for _ in range(2600):
                    line = fp_in.readline()
            else:
                fp_out.write(line)
                line = fp_in.readline()
                fp_out.write(line)
        elif line.startswith("ITEM: NUMBER OF ATOMS"):
            fp_out.write(line)
            line = fp_in.readline()
            fp_out.write(line)
            numAtoms = int(line)
            atoms = {}
        elif line.startswith("ITEM: ATOMS"):
            fp_out.write(f"ITEM: ATOMS id type xu yu zu\n")
            for _ in range(numAtoms):
                item = fp_in.readline().split()
                atoms[int(item[0])] = (item[0], item[1], item[5], item[6], item[7])
            for i in range(1, numAtoms + 1):
                temp = atoms[i]
                fp_out.write(f"{temp[0]} {temp[1]} {temp[2]} {temp[3]} {temp[4]}\n")
        else:
            fp_out.write(line)

fp_out.close()

print("Done! >:D")
