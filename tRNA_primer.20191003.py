import numpy as np

# Define the maximum number of mismatch tolerated when judging whether the RNA sequence is a tRNA.
Max_Mismatch = 2

def complement(Base): # This function return the complementary base of the input base (for RNA).
    baseset = {
        "A": "U",
        "U": "A",
        "C": "G",
        "G": "C"}
    for item in baseset:
        if Base.upper() == item:
            return baseset[item]

def complement_dna(Base): # This function return the complementary base of the input base (for DNA).
    baseset = {
        "A": "T",
        "U": "A",
        "C": "G",
        "G": "C"}
    for item in baseset:
        if Base.upper() == item:
            return baseset[item]

def complementSeq(seq): # This function return the complementary sequence of the input sequence (for RNA).
    comSeq = ""
    for i in  range(len(seq)):
        comSeq += str(complement(seq[i]))
    return comSeq

def complementSeq_dna(seq): # This function return the complementary sequence of the input sequence (for DNA).
    comSeq = ""
    for i in  range(len(seq)):
        comSeq += str(complement_dna(seq[i]))
    return comSeq

def isComplement(Base1, Base2, drawMode="none", alignMode="none"):  # This function judges whether two input bases are complementary.
    if Base1.upper() == complement(Base2):
        if alignMode.lower() == "vertical":
            return 3
        else:
            return 1
    elif (Base1.upper() == "G" and Base2.upper() == "U") or \
            (Base1.upper() == "U" and Base2.upper() == "G"):
        if drawMode.lower() == "draw":
            return 2
        else:
            return 1
    else:
        return 0

def linkerDtm(seq1, seq2, mode="none"): # This function determines the linker between two sequences.
    linkSymbol = " -+|"

    L1 = len(seq1)
    L2 = len(seq2)
    if L1 != L2:
        return 0
    else:
        L = L1
    linkerSeq = ""
    for i in range(L):
        Base1 = seq1[i]
        Base2 = seq2[L-1-i]
        linkerSeq = linkerSeq + linkSymbol[isComplement(Base1, Base2, drawMode="draw", alignMode=mode)]
    return linkerSeq

def findRegion(seq, L): # This function judges whether there is a tRNA sequence contained in the input sequence
    # and return the basic components of tRNA (e.g. TArm, TLoop, etc.) for subsequent plotting step.
    for n0 in range(L):
        for interval in range(69,94):
            # Test for Acceptor stem formation
            AcptStemStart = n0
            AcptStemEnd = n0 + interval
            AcptStemLen = 7
            MisMatch_total = 0

            if AcptStemEnd >= L:
                continue
            mismatch1 = 0
            S = 0
            for i in range(AcptStemLen):
                Base1 = seq[AcptStemStart + i]
                Base2 = seq[AcptStemEnd - i]
                S += isComplement(Base1, Base2)
            mismatch1 += (AcptStemLen - S)

            if S < AcptStemLen - Max_Mismatch or mismatch1 > Max_Mismatch:
                continue

            # Test for T arm formation
            TArmSiteEnd = n0 + interval - AcptStemLen
            TArmSiteStart = TArmSiteEnd - 16
            TArmStemLen = 5

            S = 0
            for i in range(TArmStemLen):
                Base1 = seq[TArmSiteStart + i]
                Base2 = seq[TArmSiteEnd - i]
                S += isComplement(Base1, Base2)
            mismatch1 += (TArmStemLen - S)
            MisMatch_total = mismatch1

            if S < 4 or MisMatch_total > Max_Mismatch:
                continue

            # Test for D loop formation
            for DLoopLen in range(14,19):
                MisMatch_total = mismatch1
                DLoopStart = n0 + 9
                DLoopEnd = DLoopStart + DLoopLen

                S = 0
                if DLoopLen < 17:
                    DLoopStemLen = 3
                else:
                    DLoopStemLen = 4
                for i in range(DLoopStemLen):
                    Base1 = seq[DLoopStart + i]
                    Base2 = seq[DLoopEnd - i]
                    S += isComplement(Base1, Base2)
                mismatch2 = DLoopStemLen - S
                MisMatch_total += mismatch2

                if MisMatch_total <= Max_Mismatch:
                    # Test for Anticodon loop // Problem occurs here
                    AtCLoopStart = DLoopEnd + 2
                    AtCLoopStemLen = 5
                    for AtCLoopEnd in range(max(AtCLoopStart + 16, AcptStemEnd - 54),
                                            min(AtCLoopStart + 16, AcptStemEnd - 28)+1):
                        S = 0
                        MisMatch_total = mismatch1 + mismatch2
                        for i in range(AtCLoopStemLen):
                            Base1 = seq[AtCLoopStart + i]
                            Base2 = seq[AtCLoopEnd - i]
                            S += isComplement(Base1, Base2)
                        mismatch3 = AtCLoopStemLen - S
                        MisMatch_total += mismatch3

                        if MisMatch_total <= Max_Mismatch:
                            # Print the result according to mode set by the user
                            print("%d..%d" % (AcptStemStart+1, AcptStemEnd+1+1), end='\t')
                            print("mismatch=%d" % MisMatch_total, end='\t')

                            GroupNum = [AcptStemStart, AcptStemLen, 2, 4, DLoopLen-2*4+1, 4, 1, AtCLoopStemLen,
                                        AtCLoopEnd-AtCLoopStart-2*AtCLoopStemLen+1, AtCLoopStemLen,
                                        TArmSiteStart-AtCLoopEnd-1, TArmStemLen, 16-2*TArmStemLen+1,
                                        TArmStemLen, AcptStemLen, len(seq)-AcptStemLen]
                            GroupNum = np.cumsum(GroupNum)
                            Output = {
                                "AcceptArm": "",
                                "Interval1": "",
                                "DArm":"",
                                "DLoop":"",
                                "DArmComp":"",
                                "Interval2":"",
                                "AnticodonArm":"",
                                "AnticodonLoop":"",
                                "AnticodonArmComp":"",
                                "VariableLoop":"",
                                "TArm":"",
                                "TLoop":"",
                                "TArmComp":"",
                                "AcceptArmComp":"",
                                "Redundance":""
                            }
                            OutputName = list(Output.keys())
                            for i in range(AcptStemStart, len(seq)):
                                for j in range(len(GroupNum)-1):
                                    if i in range(GroupNum[j], GroupNum[j+1]):
                                        Output[OutputName[j]] = Output[OutputName[j]] + seq[i]

                            return Output
    print("Whoops! tRNA judge thinks it's not a tRNA. Please check...")

def draw_tRNA(ForDraw, mode):   # This function draws a tRNA.
    if mode == 0:
        print("Template version...", end='\n\n')
    if mode == 1:
        print("Designed version...", end='\n\n')

    # ----------------------------- Revision start -----------------------------
    print("                        %s-5" % ForDraw["Redundance"])
    print("                 3-%c %c %c" % (ForDraw["AcceptArm"][0],
                                        linker["AcceptorArmLink"][0],
                                        ForDraw["AcceptArmComp"][len(ForDraw["AcceptArm"]) - 1]))
    # ----------------------------- Revision end -----------------------------

    for i in range(1,len(ForDraw["AcceptArm"])-1):
        print("                   %c %c %c" % (ForDraw["AcceptArm"][i],
                                             linker["AcceptorArmLink"][i],
                                             ForDraw["AcceptArmComp"][len(ForDraw["AcceptArm"])-1-i]))
    print("                   %c %c %c          %c %c" % (ForDraw["AcceptArm"][len(ForDraw["AcceptArm"])-1],
                                                        linker["AcceptorArmLink"][len(ForDraw["AcceptArm"])-1],
                                                        ForDraw["AcceptArmComp"][0],
                                                        ForDraw["TLoop"][6],
                                                        ForDraw["TLoop"][5]))
    print("                  %c     %c %c %c %c %c     %c" % (ForDraw["Interval1"][0],
                                                            ForDraw["TArmComp"][4],
                                                            ForDraw["TArmComp"][3],
                                                            ForDraw["TArmComp"][2],
                                                            ForDraw["TArmComp"][1],
                                                            ForDraw["TArmComp"][0],
                                                            ForDraw["TLoop"][4]))
    print("      %c %c        %c      %c %c %c %c %c     %c" % (ForDraw["DLoop"][1],
                                                               ForDraw["DLoop"][0],
                                                               ForDraw["Interval1"][1],
                                                               linker["TArmLink"][0],
                                                               linker["TArmLink"][1],
                                                               linker["TArmLink"][2],
                                                               linker["TArmLink"][3],
                                                               linker["TArmLink"][4],
                                                               ForDraw["TLoop"][3]))
    print("    %c     %c %c %c %c       %c %c %c %c %c     %c" % (ForDraw["DLoop"][2],
                                                                   ForDraw["DArm"][3],
                                                                ForDraw["DArm"][2],
                                                                ForDraw["DArm"][1],
                                                                ForDraw["DArm"][0],
                                                                ForDraw["TArm"][0],
                                                                ForDraw["TArm"][1],
                                                                ForDraw["TArm"][2],
                                                                ForDraw["TArm"][3],
                                                                ForDraw["TArm"][4],
                                                                ForDraw["TLoop"][2]))

    # ----------------------------- Revision start -----------------------------
    print("    %c     %c %c %c %c        %c        %c %c" % (ForDraw["DLoop"][3],
                                                           linker["DArmLink"][0],
                                                           linker["DArmLink"][1],
                                                           linker["DArmLink"][2],
                                                           linker["DArmLink"][3],
                                                           ForDraw["VariableLoop"][4],
                                                           ForDraw["TLoop"][1],
                                                           ForDraw["TLoop"][0]))
    if Type_RT == "MMLV" or Type_RT == "RSV":
        print("      %c %c %c %c %c %c         %c" % (ForDraw["DLoop"][4],
                                                        ForDraw["DLoop"][6],
                                                        ForDraw["DArmComp"][0],
                                                        ForDraw["DArmComp"][1],
                                                        ForDraw["DArmComp"][2],
                                                        ForDraw["DArmComp"][3],
                                                        ForDraw["VariableLoop"][3]))
        print("                 %c     %c %c %c" % (ForDraw["Interval2"],
                                                    ForDraw["VariableLoop"][0],
                                                    ForDraw["VariableLoop"][1],
                                                    ForDraw["VariableLoop"][2]))
    elif Type_RT == "HIV1":
        print("    %c     %c %c %c %c         %c  " % (ForDraw["DLoop"][4],
                                                       ForDraw["DArmComp"][0],
                                                       ForDraw["DArmComp"][1],
                                                       ForDraw["DArmComp"][2],
                                                       ForDraw["DArmComp"][3],
                                                       ForDraw["VariableLoop"][4]))
        print("      %c %c        %c     %c %c %c" % (ForDraw["DLoop"][6],
                                                      ForDraw["DLoop"][7],
                                                      ForDraw["Interval2"],
                                                      ForDraw["VariableLoop"][0],
                                                      ForDraw["VariableLoop"][1],
                                                      ForDraw["VariableLoop"][2],))
    # ----------------------------- Revision end -----------------------------

    for i in range(len(ForDraw["AnticodonArm"])):
        print("                  %c %c %c" % (ForDraw["AnticodonArm"][i],
                                              linker["AnticodonArmLink"][i],
                                              ForDraw["AnticodonArmComp"][len(ForDraw["AnticodonArm"])-1-i]))
    print("                %c       %c " % (ForDraw["AnticodonLoop"][0],
                                            ForDraw["AnticodonLoop"][6]))
    print("                %c       %c " % (ForDraw["AnticodonLoop"][1],
                                            ForDraw["AnticodonLoop"][5]))
    print("                  %c %c %c" % (ForDraw["AnticodonLoop"][2],
                                          ForDraw["AnticodonLoop"][3],
                                          ForDraw["AnticodonLoop"][4]))

## Main

# Select the type of reverse transcriptase and the corresponding template tRNA is defined.
Type_RT = str(input("Please input the type of reverse transcriptase that you want to use (select from MMLV/HIV1/RSV): "))

if Type_RT == "MMLV":
    temp_seq = "GGCUCGUUGGUCUAGGGGUAUGAUUCUCGCUUAGGGUGCGAGAGGUCCCGGGUUCAAAUCCCGGACGAGCCCCCA"
    PBS_len = 18
elif Type_RT == "HIV1":
    temp_seq = "GCCCGGCUAGCUCAGUCGGUAGAGCAUCAGACUUUUAAUCUGAGGGUCCAGGGTUCAAGUCCCUGUUCGGGCGCCA"
    PBS_len = 18
elif Type_RT == "RSV":
    temp_seq = "GACCUCGUGGCGCAACGGUAGCGCGUCUGACUCCAGAUCAGAAGGCUGCGUGUUCGAAUCACGUCGGGGUCACCA"
    PBS_len = 17
else:
    exit("Invalid reverse transcriptase type. Please check.")

# Input the sequence of gene of interest.
Mut_seq = str(input("Please input a sequence that you want to mutate: "))
# Mut_seq = "CGATGATCGATCAGCTAGCTAGTCGATCGATCGATC"

if len(Mut_seq) < PBS_len:
    exit("The length of the sequence should be larger than or equal to 18nt for MMLV RT and HIV-1 RT and should be larger than or equal to 17nt for RSV RT. Please check your input.")

for nt in Mut_seq:
    if nt != "A" and nt != "T" and nt != "C" and nt != "G":
        exit("A DNA sequence should only contain character A/T/C/G. Please check your input.")

# Derive the PBS sequence as well as the complementary sequence that should be used to replace several 3'-terminal nucleotides of template tRNA.
PBS = Mut_seq[len(Mut_seq)-1:len(Mut_seq)-PBS_len-1:-1]
Primerseq = ""
for i in range(len(PBS)):
    Primerseq += str(complement(PBS.replace("T","U")[i]))

# Draw template tRNA (tRNA_Pro for MMLV RT, tRNA_Lys for HIV-1 RT, tRNA_Trp for RSV RT).
L_temp = len(temp_seq)
ForDraw = findRegion((temp_seq.upper()).replace("T","U"), L_temp)
if ForDraw != None:

    if len(ForDraw["DLoop"]) < 9:
        for i in range(9 - len(ForDraw["DLoop"])):
            ForDraw["DLoop"] = ForDraw["DLoop"] + " "

    linker = {
        "AcceptorArmLink":"",
        "DArmLink":"",
        "AnticodonArmLink":"",
        "TArmLink":""
    }

    linker["AcceptorArmLink"] = linkerDtm(ForDraw["AcceptArm"], ForDraw["AcceptArmComp"])
    linker["DArmLink"] = linkerDtm(ForDraw["DArmComp"],ForDraw["DArm"], mode='vertical')
    linker["AnticodonArmLink"] = linkerDtm(ForDraw["AnticodonArm"], ForDraw["AnticodonArmComp"])
    linker["TArmLink"] = linkerDtm(ForDraw["TArm"], ForDraw["TArmComp"], mode='vertical')

    draw_tRNA(ForDraw, mode=0)

# Replace several 3'-terminal nucleotides of template tRNA with primer sequence that should be complementary to PBS.
ForDraw["Redundance"] = Primerseq[len(ForDraw["Redundance"])-1::-1]
ForDraw["AcceptArmComp"] = Primerseq[len(ForDraw["Redundance"])+len(ForDraw["AcceptArmComp"])-1:len(ForDraw["Redundance"])-1:-1]
ForDraw["AcceptArm"] = complementSeq(ForDraw["AcceptArmComp"])[::-1]
ForDraw["TArmComp"] = Primerseq[len(ForDraw["Redundance"])+len(ForDraw["AcceptArmComp"]):len(ForDraw["Redundance"])+len(ForDraw["AcceptArmComp"])+5]
ForDraw["TArmComp"] = ForDraw["TArmComp"][::-1]
ForDraw["TArm"] = complementSeq(ForDraw["TArmComp"])[::-1]
remain = Primerseq[len(ForDraw["TArmComp"])+len(ForDraw["Redundance"])+len(ForDraw["AcceptArmComp"])::][::-1]
ForDraw["TLoop"] = ForDraw["TLoop"][:len(ForDraw["TLoop"])-len(remain)] + remain

# Draw tRNA primer
if ForDraw != None:

    if len(ForDraw["DLoop"]) < 9:
        for i in range(9 - len(ForDraw["DLoop"])):
            ForDraw["DLoop"] = ForDraw["DLoop"] + " "

    linker = {
        "AcceptorArmLink":"",
        "DArmLink":"",
        "AnticodonArmLink":"",
        "TArmLink":""
    }

    linker["AcceptorArmLink"] = linkerDtm(ForDraw["AcceptArm"], ForDraw["AcceptArmComp"])
    linker["DArmLink"] = linkerDtm(ForDraw["DArmComp"],ForDraw["DArm"], mode='vertical')
    linker["AnticodonArmLink"] = linkerDtm(ForDraw["AnticodonArm"], ForDraw["AnticodonArmComp"])
    linker["TArmLink"] = linkerDtm(ForDraw["TArm"], ForDraw["TArmComp"], mode='vertical')

    # Draw template tRNA
    draw_tRNA(ForDraw, mode=1)

# Print the DNA sequence encoding the tRNA primer
print("Your tRNA primer: %s%s%s%s%s%s%s%s%s%s%s%s%s" % (complementSeq_dna(ForDraw["AcceptArm"]),
                                                        complementSeq_dna(ForDraw["Interval1"]),
                                                        complementSeq_dna(ForDraw["DArm"]),
                                                        complementSeq_dna(ForDraw["DArmComp"]),
                                                        complementSeq_dna(ForDraw["Interval2"]),
                                                        complementSeq_dna(ForDraw["AnticodonArm"]),
                                                        complementSeq_dna(ForDraw["AnticodonLoop"]),
                                                        complementSeq_dna(ForDraw["VariableLoop"]),
                                                        complementSeq_dna(ForDraw["TArm"]),
                                                        complementSeq_dna(ForDraw["TLoop"]),
                                                        complementSeq_dna(ForDraw["TArmComp"]),
                                                        complementSeq_dna(ForDraw["AcceptArmComp"]),
                                                        complementSeq_dna(ForDraw["Redundance"])), end="\n\n")

print("End...")
