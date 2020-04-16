# Julie Shay
# April 1, 2020
# This script will convert a CIPARS table into an NCBI antibiogram submission table.
import pandas as pd
import argparse

ncbi_columns = ["sample_name/biosample_accession", "antibiotic", "resistance_phenotype", "measurement_sign", "measurement", "measurement_units",  "laboratory_typing_method", "laboratory_typing_platform", "vendor", "laboratory_typing_method_version_or_reagent", "testing_standard"]
cipars_constants = {"laboratory_typing_method": "MIC", "laboratory_typing_platform": "Sensititre", "measurement_units": "mg/L"}
# process a file in the standard CIPARS format
def processCIPARS(infile):
    xls = pd.ExcelFile(infile)
    # grab information about the panel/break points
    try:
        # never mind about grabbing the panel version--there's no place to put it in the antibiogram table
        # panel = pd.read_excel(xls, "Break Points", skiprows=4, nrows=0).columns[0][12:]
        paneldf = pd.read_excel(xls, "Break Points", skiprows=5, skipfooter=8, usecols=["Antimicrobial", "Susceptible", "Resistant"])
        paneldf["testing_standard"] = paneldf["Antimicrobial"].str[-1:]
        paneldf["Antimicrobial"] = paneldf["Antimicrobial"].str[0:-2]
        paneldf = paneldf.set_index("Antimicrobial")
        paneldf["Susceptible"] = paneldf["Susceptible"].str.replace("≤", "").astype(float)
        paneldf["Resistant"] = paneldf["Resistant"].str.replace("≥", "").astype(float)
    except:
        exit("Input file did not contain a Break Points tab in the expected CIPARS format.")

    # now try to figure out what the guidelines are...it just free text, so this probably won't be robust...
    try:
        guidelines = pd.read_excel(xls, "Break Points", skiprows=(6 + len(paneldf)), usecols=[0], nrows=len(paneldf["testing_standard"].unique()), names=["value"])
        guidelines["key"] = guidelines["value"].str[0]
        guidelines["value"] = guidelines["value"].str[1:]
        guidelines = guidelines.set_index("key").to_dict()["value"]
        paneldf["testing_standard"] = paneldf["testing_standard"].map(guidelines)
    except:
        print("Couldn't figure out the testing_standard used.")
    
    # grab actual AMR data
    outdf = pd.DataFrame(columns=ncbi_columns)
    for sheet in xls.sheet_names[2:]:
        datadf = pd.read_excel(xls, sheet)
        for _, row in datadf.iterrows():
            outdf = outdf.append(CIPARSline(row, paneldf))
    for key in cipars_constants:
        outdf[key] = cipars_constants[key]
    return outdf


# process information about a single sample from a file in the standard CIPARS format
def CIPARSline(inline, antibiotics):
    amrtable = pd.merge(antibiotics, inline.to_frame().reset_index(), left_on="Antimicrobial", right_on="index", how="inner")
    mcol = amrtable.columns[-1]
    amrtable["measurement_sign"], amrtable["measurement"] = splitmeasurement(amrtable[mcol])
    amrtable["resistance_phenotype"] = amrtable.apply(define_susceptible, axis=1)
    amrtable = amrtable.drop([mcol, "Susceptible", "Resistant"], axis=1).rename(columns={"Antimicrobial": "antibiotic", "index": "antibiotic"})
    amrtable["sample_name/biosample_accession"] = inline["Sample ID"]
    return amrtable


# probably far from the most efficient way to split up a measurement value and a measurement sign
def splitmeasurement(inseries):
    try:
        df = inseries.str[::-1].str.split(" ", expand=True)
        return df[1].str[::-1], df[0].str[::-1].astype(float)
    except:
        return None, inseries.astype(float)


# just figure out whether the "measurement" is considered resistant/susceptiblee
def define_susceptible(i):
    try:
        i["measurement"] = float(i["measurement"])
        if i["measurement"] <= i["Susceptible"]:
            return "susceptible"
        elif i["measurement"] >= i["Resistant"]:
            return "resistant"
        else:
            return "intermediate"
    except:
        return "not defined"


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', dest="infile", type=str, required=True, help="""Input Excel file""")
    parser.add_argument('-o', dest="outfile", default="ncbi_amr.tsv", help="""Output TSV file""")
    args = parser.parse_args()
    outdf = processCIPARS(args.infile)
    outdf.to_csv(args.outfile, sep="\t", index=False)


if __name__ == "__main__":
    main()

