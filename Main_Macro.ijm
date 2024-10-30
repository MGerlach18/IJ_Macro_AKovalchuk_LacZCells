/*
 * Macro to processWF Images of DAPI/LacZ staines cells from ApoTome
 * Automatic classification and quantification of nuclei
 */

//getting input parameters
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ Float (label = "Advanced: Distance Nuclei-LacZ", style = "slider", min=0.1, max=4, stepSize=0.1, value=1) distance


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

//function to open and process Images  
function processFile(input, output, file) {
//BioFormat extension to find the number of positions in a file
Ext.setId(input + File.separator + list[i]);
Ext.getSeriesCount(seriesCount);

for (series = 1; series <= seriesCount; series++) {
run("Bio-Formats Importer", "open=[" + input + File.separator + list[i] + "] color_mode=Default view=Hyperstack stack_order=XYCZT series_"+series");
title=getTitle();
run("Split Channels");

/*
 * Segment LacZ Signal
 * create psrticles with presettings
 * create mask from selection
 * Distance map
 */
 
 /*
  * Segment DAPI channel
  * dectect particles to ROI manager
  * Lopp ROI --> Distance to LacZ signal
  * counter for nuclei, classes
  * report to result table
  * Save ROIS
  */
roiManager("Save", output + "\\ROIS\\" + title + ".zip");
roiManager("Deselect");
roiManager("Delete");

selectWindow("Summary_Total");
for (o = 0; o < count; o++) {
Table.set("Name", o+a, title);
Table.set("ROI #", o+a, o);
Table.set("Area [µm²]", o+a, getResult("Area", o));
}

}






close("Results");
roiManager("Deselect");
roiManager("Delete");
close("*");
}
