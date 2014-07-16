import Migration
import os

# define the paths where the files are located
configpath="/Volumes/TimeMachineBackups/Mariana_all_confocal/Mariana_migration_output/sexconfig/"

# inpath="/Volumes/TimeMachineBackups/Mariana_all_confocal/Mariana_migration_adam_copy/day_1/no_tgfb/gel_1/ZSeries-08012013-0901-006/dapi/"
# outpath="/Volumes/TimeMachineBackups/Mariana_all_confocal/Mariana_migration_output/day_1/no_tgfb/gel_1/trial_4/"
def run(inpath,outpath,configpath=configpath):
    # assert paths exist
    # if not os.path.exists(inpath):
    #     inpath = inpath.replace('/dapi','/') # if dapi doesnt exist, try the parent directory
    #     
    # initialize the stack
    imgstack = Migration.LSMStack(image_directory=inpath,output_directory=outpath,config_directory=configpath)

    # convert the images to .fits files
    imgstack.PrepareImages()

    # run sextractor on the images (set checkimages=True for additional output images)
    imgstack.FindDetections(checkimages=True)

    # read in the sextractor output into a table
    imgstack.ReadCatalogs()

    # run the object/association assignment algorithm
    imgstack.FindObjects()

    # save the database
    imgstack.Save()

    # optional: plot count histogram
    imgstack.PlotCountHist()

    # optional: export to an excel file 
    imgstack.ExportToExcel()

    # optional: prepare verification images (takes several minutes)
    # imgstack.PrepareVerificationImages()

    # import glob
    # pklpath = glob.glob(outpath+'*.pkl')[0]
    # imgstack = Migration.loadPickle(pklpath)

    # optional: create verification images
    # imgstack.MakeVerificationImages()


