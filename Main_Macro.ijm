
//Splitting channels
run("Split Channels");
close("C3-" + title);

//Processing LacZ Fluorescence channel to detect LacZ patches
selectWindow("C1-" + title);
run("Subtract Background...", "rolling=100");
run("Gaussian Blur...", "sigma=3");
setAutoThreshold("Intermodes dark");
run("Convert to Mask");
run("Analyze Particles...", "size=20-Infinity show=Masks include");
run("Invert");
run("Chamfer Distance Map", "distances=[Quasi-Euclidean (1,1.41)] output=[32 bits] normalize");
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
}
//cleaning up workspace for next series

}

