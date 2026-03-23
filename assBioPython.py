from functionsBio import bioFun
import functionsRun

opt,args = functionsRun.fillingOptsArgs()
mybio = bioFun(args[1])
functionsRun.runArgs(mybio,args,opt)

# python assBioPython dna AGCTGACTGACTACGTCGAGTCGTACGCA
