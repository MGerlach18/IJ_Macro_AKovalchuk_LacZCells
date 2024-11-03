/*
 * Macro to processWF Images of DAPI/LacZ staines cells from ApoTome
 * Automatic classification and quantification of nuclei
 * Requires: 
 * ImageJ version 1.54f
 * MorphoLibJ
 * BioFormats Extension
 * CC-BY 4.0 by Michael Gerlach, TU Dresden
 */

//getting input parameters
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ Float (label = "Advanced: Distance Nuclei-LacZ [µm]", style = "slider", min=0, max=5, stepSize=0.1, value=1) distance


//Preparing Stage
run("Bio-Formats Macro Extensions");
setOption("BlackBackground", true);
run("Set Measurements...", "min redirect=None decimal=0")
print("\\Clear");
run("Clear Results");
close("*");
setBatchMode(true);

if (roiManager("count")>0) {
roiManager("Deselect");
roiManager("Delete");
}

File.makeDirectory(output + "\\ROIS");
File.makeDirectory(output + "\\Results");

//Variable for internal counting
var line=0;

//variables for Metadata
var plate 
var well
var position

//measuring the folder and creating results table
Table.create("Summary_Total");
//function to create cumulated stack from all files in the input folder
processFolder(input);

//Processing the cumulated stack
selectWindow("ExperimentStack");
stack.getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);

//Splitting channels
run("Split Channels");
//close("C3-" + title);

//Processing LacZ Fluorescence channel to detect LacZ patches
selectWindow("C1-" + title);
run("Subtract Background...", "rolling=100");
run("Gaussian Blur...", "sigma=3");
setAutoThreshold("Intermodes dark stack");
run("Convert to Mask", "method=Intermodes background=Dark calculate black");
run("Analyze Particles...", "size=20-Infinity show=Masks include stack");
run("Invert", "stack");
run("Distance Map", "stack");
rename("DistanceMap");

//Processing Hoechst Fluorescence channel to detect nuclei 
selectWindow("C2-" + title);
run("Gaussian Blur...", "sigma=3");
setAutoThreshold("Default dark");
run("Convert to Mask");
run("Analyze Particles...", "  show=Nothing exclude include add");
close();
 
//Measuring distance of each Nucleus according to detected LacZ patches and assigning it to groups 1 (positive) or 2 (negative) --> internal counting
selectWindow("DistanceMap");
roiManager("deselect");
roiManager("measure");
positive=0;
negative=0;
n = roiManager('count');
for (a = 0; a < n; a++) {
    roiManager('select', a);
    m=getResult("Min", a);
    if (m <= (distance/pixelWidth)) {
    	positive=positive+1;	
    	RoiManager.setGroup(1);
    }
    else{
    	negative=negative+1;
    	RoiManager.setGroup(2);	
    }
  
}
  run("Clear Results");
//creating entry to results table - one line per position

selectWindow("Summary_Total");
Table.set("Plate Name", line, list[i]);
Table.set("Plate well", line, Well);
Table.set("Position number", line, number);
Table.set("Positive nuclei", line, positive);
Table.set("Negative nuclei", line, negative);
Table.set("% Pos", line, positive/(positive+negative)*100);
Table.update;
line=line+1;

roiManager("save", output + File.separator + "ROIS" + File.separator + list[i]+"_"+SeriesName+".zip");
roiManager("Deselect");
roiManager("Delete");
close("*");



for (u = 0; i < slices; u++) {


//save Results as .csv-file and clean up
selectWindow("Summary_Total");
saveAs("Results", output + "\\Results\\Results_" + distance + "µm.csv");
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
			processFile(input, output, list[i]);
	}
}

//function to open and process Images  
function processFile(input, output, file) {
//BioFormat extension to find the number of positions in a file
Ext.setId(input + File.separator + list[i]);
Ext.getSeriesCount(seriesCount);


//Process all postions in a file
for (series = 1; series <= seriesCount; series++) {
//Detect Well Postion and Position number in .CZI file metadata
SeriesName=IJ.pad(series, 2);
plate[line]=SeriesName;

Ext.getMetadataValue("Information|Image|S|Scene|ArrayName #" + SeriesName, Well);
well[line]=Well;
Ext.getMetadataValue("Information|Image|S|Scene|Name #" + SeriesName, number);
position[line]=number;

run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+series);
title=getTitle();

if {line==0}{
	rename("ExperimentStack");
	line=line+1;
} else {
	run("Concatenate...", " title=ExperimentStack open image1=ExperimentStack image2=" + title);
	line=line+1;
}

}
}
