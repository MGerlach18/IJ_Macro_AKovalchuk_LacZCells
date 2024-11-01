/*
 * Macro to processWF Images of DAPI/LacZ staines cells from ApoTome
 * Automatic classification and quantification of nuclei
 * Requires: 
 * ImageJ version
 * MorphoLibJ
 * BioFormats Extension
 * CC-BY 4.0 by Michael Gerlach, TU Dresden
 */

//getting input parameters
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ Float (label = "Advanced: Distance Nuclei-LacZ", style = "slider", min=0.1, max=4, stepSize=0.1, value=1) distance


//Preparing Stage
run("Bio-Formats Macro Extensions");
setOption("BlackBackground", true);
run("Set Measurements...", "min redirect=None decimal=0")
print("\\Clear");
run("Clear Results");
close("*");
setBatchMode(false);

if (roiManager("count")>0) {
roiManager("Deselect");
roiManager("Delete");
}

File.makeDirectory(output + "\\ROIS");
File.makeDirectory(output + "\\Results");

//measuring the folder and creating results table
Table.create("Summary_Total");
processFolder(input);

//save Results as .csv-file and clean up
selectWindow("Summary_Total");
saveAs("Results", output + "\\Results\\Results.csv");
close("*");
print("Batch processing completed");

//End of Macro




//Definition of functions below

// functions to scan folders/subfolders/files to find files with correct suffix and do the (pre-)processing

function processFolder(input) {
	list = getFileList(input);
	for (i = 0; i < list.length; i++) {
		if(File.isDirectory(input + File.separator + list[i]))
			processFolder(input + File.separator + list[i]);
		if(endsWith(list[i], ".czi"))
			measureFile(input, output, list[i]);
	}
}

//function to open and process Images  
function processFile(input, output, file) {
//BioFormat extension to find the number of positions in a file
Ext.setId(input + File.separator + list[i]);
Ext.getSeriesCount(seriesCount);

for (series = 1; series <= seriesCount; series++) {
run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+series");
title=getTitle();
run("Split Channels");
selectWindow("C1-" + title);
run("Subtract Background...", "rolling=100");
run("Gaussian Blur...", "sigma=3");
setAutoThreshold("Intermodes dark");
run("Convert to Mask");
run("Analyze Particles...", "size=20-Infinity show=Masks include");
run("Invert");
run("Chamfer Distance Map", "distances=[Quasi-Euclidean (1,1.41)] output=[32 bits] normalize");
ID=getImageID();
 
 
selectWindow("C2-" + title);
run("Gaussian Blur...", "sigma=3");
setAutoThreshold("Default dark");
run("Convert to Mask");
run("Analyze Particles...", "  show=Nothing exclude include add");
close();
 
selectImage(ID);
n = roiManager('count');
roiManager("measure");
positive=0;
negative=0;
for (a = 0; a < n; a++) {
    roiManager('select', a);
    m=getResult("Min", a);
    if (m < (distance*pixelWidth)) {
    	positive=positive+1;	
    	RoiManager.setGroup(1);
    }
    else{
    	negative=negative+1;
    	RoiManager.setGroup(2);	
    }
    
run("Clear Results");
selectWindow("Summary_Total");
Table.set("Name", i+series, title);
Table.set("Positive nuclei", i+series, positive);
Table.set("Negative nuclei", i+series, negative);
Table.set("% Pos", i+series, positive/(positive+negative)*100);
}
roiManager("Deselect");
roiManager("Delete");
close("*");
}
}
