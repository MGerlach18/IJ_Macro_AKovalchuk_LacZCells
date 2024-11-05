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
setBatchMode(false);

if (roiManager("count")>0) {
roiManager("Deselect");
roiManager("Delete");
}

File.makeDirectory(output + "\\ROIS");
File.makeDirectory(output + "\\Results");
File.makeDirectory(output + "\\Results\\Controls");

//Variable for internal counting
var line=0;

//variables for Metadata
var plate = newArray();
var well = newArray();
var position = newArray();

//creating results table
Table.create("Summary_Total");

//create cumulated stack from all files in the input folder
processFolder(input);

//Processing the cumulated stack
selectWindow("ExperimentStack");
getDimensions(width, height, channels, slices, frames);
getPixelSize(unit, pixelWidth, pixelHeight);

//Splitting channels
run("Split Channels");
//close("C3-" + title);

//Processing LacZ Fluorescence channel to detect LacZ patches
selectWindow("C1-ExperimentStack");
run("Subtract Background...", "rolling=100 stack");
run("Gaussian Blur...", "sigma=3 stack");
setAutoThreshold("Moments dark no-reset stack");
run("Convert to Mask", "method=Moments background=Dark black list");
run("Analyze Particles...", "size=20-Infinity show=Masks include stack");
selectWindow("Mask of C1-ExperimentStack");
//run("Invert", "stack");
run("Distance Map", "stack");
rename("DistanceMap");

//Processing Hoechst Fluorescence channel to detect nuclei 
selectWindow("C2-ExperimentStack");
run("Gaussian Blur...", "sigma=3 stack");
setAutoThreshold("Moments dark no-reset stack");
run("Convert to Mask", "method=Default background=Dark black");
run("Watershed", "stack");
run("Analyze Particles...", "  show=Nothing exclude include add stack");
RoiManager.associateROIsWithSlices(true);
close();
 
//Measuring distance of each Nucleus according to detected LacZ patches and assigning it to groups 1 (positive) or 2 (negative) --> internal counting
for (u = 1; u <= line; u++) {
	positive=0;
	negative=0;
	n = roiManager('count');
	for (x = 0; x < n; x++) {
   	 roiManager('select', x);
   	 RoiName=Roi.getName;
   	 RoiSlice=split(RoiName, "-");
   	 if (RoiSlice[0]==IJ.pad(u, 4)) {
   	 	run("Measure");
   	 	m=getResult("Min", 0);
   	 	if (m <= (distance/pixelWidth)) {
    	positive=positive+1;	
    	RoiManager.setGroup(1);
   		 }
   		 else{
    			negative=negative+1;
    			RoiManager.setGroup(2);	
    		}
   	 }
	}
   	selectWindow("Summary_Total");
   	Table.set("Index", u-1, u);
	Table.set("Plate Name", u-1, plate[u-1]);
	Table.set("Plate well", u-1, well[u-1]);
	Table.set("Position number", u-1, position[u-1]);
	Table.set("Positive nuclei", u-1, positive);
	Table.set("Negative nuclei", u-1, negative);
	Table.set("% Pos", u-1, positive/(positive+negative)*100);
	Table.update;
   	}

//SavingRoiSet
roiManager("save", output + File.separator + "ROIS" + File.separator + "ROI_Set.zip");

//Creating Control Previews
selectWindow("C2-ExperimentStack");
run("Enhance Contrast", "saturated=0.35");
roiManager("Show All");
run("Flatten", "stack");
run("Image Sequence... ", output + File.separator + "Controls\\ format=JPEG name=Control_"+distance+" use");

roiManager("Deselect");
roiManager("Delete");
close("*");

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

if (line==0){
	rename("ExperimentStack");
	line=line+1;
} else {
	run("Concatenate...", " title=ExperimentStack open image1=ExperimentStack image2=[" + title + "]");
	line=line+1;
}

}
}
