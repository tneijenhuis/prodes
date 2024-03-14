import argparse
import json
from prodes import data

def parse_arguments():
    """parses the arguments"""

    parser = argparse.ArgumentParser(description="reads pka output files and returns a dictionary or json formatted file")
    parser.add_argument("pka_file", help="Path to the pka file to parse")
    parser.add_argument("source", help="software used to generate the pka file")
    parser.add_argument("-o", "--output", help="Path to the output file ")
    arg = parser.parse_args()

    return arg.pka_file, arg.source, arg.output

def write_json(dictionary, outputfile):
    """writes dictionary as a json"""

    written_json = json.dumps(dictionary, indent=4)

    with open(outputfile, "w") as f:
        f.write(written_json)

        
def convert_hpp(pka_file):
    """converts H++ format to dictionary"""

    pkas = {}
    with open (pka_file) as f:
        for i, line in enumerate(f):
            
            if i == 0:
                pass
            
            else:
                if line[0:4] == "Site":
                    break

                else:
                    splitted_line = line.split()
                    res_numb = splitted_line[0].split("-")[1]
                    identifier = splitted_line[0].split("-")[0]
                    res_pka = splitted_line[2]
                    if ">" in res_pka:
                        res_pka = 14.0

                    elif "<" in res_pka:
                        res_pka = 0

                    res_pka = float(res_pka)

                    if identifier[0:2] == "NT":
                        identifier = "N+"

                    elif identifier[0:2] == "CT":
                        identifier = "C-"
                    
                    else:
                        identifier = identifier[0:3].upper()

                    if res_numb in pkas:
                        pkas[res_numb].append({identifier: res_pka})
                    else:
                        pkas[res_numb] = [{identifier: res_pka}]
    return pkas

def convert_propka(pka_file):
    """converts PROPKA format to dictionary"""

    known_residues = data.all_residues()
    summary = False
    line_numb = 1
    pkas = {}
    with open(pka_file) as f:
        for line in f:

            if "----" in line:
                summary = False

            if summary:
                if line_numb == 1:
                    pass
                else:

                    identifier = line.strip()[0:4].strip()
                    if identifier in known_residues or identifier in ["N+", "C-"]:
                        res_numb = int(line.strip()[4:7].strip())
                        res_pka = float(line.strip()[12:18].strip())

                        if res_numb in pkas:
                            pkas[res_numb].append({identifier: res_pka})
                        else:
                            pkas[res_numb] = [{identifier: res_pka}]
                line_numb += 1

            if "SUMMARY OF THIS PREDICTION" in line:
                summary = True
                
    return pkas


def convert_pypka(pka_file):
    """Converts pypka output"""

    from prodes.data import residue_data
    
    reading = False
    pkas = {}
    with open(pka_file) as f:
        for line in f:
            
            if reading:
                if line[:3] == "API":
                    break
                
                split_line = line.split()

                res_numb = split_line[1]
                identifier = split_line[2]
                res_pka = split_line[3]
                if identifier != "SER" and identifier != "THR":
                    if res_pka == "Not":
                        potential_charge = residue_data(identifier)["potential_charge"]
                        if potential_charge > 0:
                            res_pka = 14
                        else:
                            res_pka = 0

                    res_pka = float(res_pka)

                    if identifier == "NTR":
                        identifier = "N+"

                    elif identifier == "CTR":
                        identifier = "C-"

                    if res_numb in pkas:
                        pkas[res_numb].append({identifier: res_pka})
                    else:
                        pkas[res_numb] = [{identifier: res_pka}]

            else: 
                if line[:5] == "Chain":
                    reading = True
    return pkas
def main():
    pka_file, source, output_file = parse_arguments()

    convert = f"convert_{source}"
    pkas = globals()[convert](pka_file)
    print(output_file)
    if output_file:
        write_json(pkas, output_file)
    else:
        print(pkas)


if __name__ == "__main__":
    main()


