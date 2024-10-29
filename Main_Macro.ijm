/*
 * Macro to processWF Images of DAPI/LacZ staines cells from ApoTome
 * Automatic classification and quantification of nuclei
 */

//getting input parameters
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output


//Preparing Stage
run("Bio-Formats Macro Extensions");
setOption("BlackBackground", true);
run("Set Measurements...", "area min redirect=None decimal=0");
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

//function to semiautomatically annotate ROIs in images  
function processFile(input, output, file) {
//Import
run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_1");
title=getTitle();

//check whether the presegmentation for this file already exists
if (File.exists(output + "\\ROIS\\" + title  +".zip")) {
	close("*");
} else {

run("RGB Color");
run("16-bit");

run("Top Hat...", "radius=7 light");

setAutoThreshold("Triangle dark");
setOption("BlackBackground", true);
run("Convert to Mask");

run("Remove Outliers...", "radius=3 threshold=50 which=Dark");

run("Analyze Particles...", "size=30.00-Infinity circularity=0.30-1.00 exclude include add");

selectWindow(title);
close("\\Others");

//Manual check and addition of detected ROIs

roiManager("Show All");
setTool("freehand");
waitForUser("Manual correction", "Please delete invalid ROIs [DEL] and manually add undetected regions via [T]. When finished, click [OK]");

roiManager("Save", output + "\\ROIS\\" + title + ".zip");
roiManager("Deselect");
roiManager("Delete");
}
}

//function opens all ROIS, measures the crossection area and includes them into one .csv-File

function measureFile(input, output, file) {
run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_1");
title=getTitle();

roiManager("Open", output + "\\ROIS\\" + title + ".zip");
selectWindow("Summary_Total");
a=getValue("results.count");

roiManager("Measure");

selectWindow("Results");
selectWindow("Summary_Total");
for (o = 0; o < count; o++) {
Table.set("Name", o+a, title);
Table.set("ROI #", o+a, o);
Table.set("Area [µm²]", o+a, getResult("Area", o));
}

close("Results");
roiManager("Deselect");
roiManager("Delete");
close("*");
}
