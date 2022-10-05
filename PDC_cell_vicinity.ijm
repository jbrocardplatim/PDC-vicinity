//Macro used to measure number of multi-labelled cells, within 5 or 10 um of PDC cells,
//from .czi stacks of five channels located in a directory to be chosen
//@Jacques Brocard for Marlène Dreux (CIRI), 2022

//--- VARIABLES needed for final classification : 
//w/ or w/o PDC in 5 um; w/ or w/o PDC in 10 um; 7 counts for each
var tot=newArray();
var nbpdc5=newArray();
var pdc5_c1=newArray();
var pdc5_c2=newArray();
var pdc5_c3=newArray();
var pdc5_c4=newArray();
var pdc5_c5=newArray();
var pdc5_c6=newArray();
var pdc5_c7=newArray();
var nbnopdc5=newArray();
var nopdc5_c1=newArray();
var nopdc5_c2=newArray();
var nopdc5_c3=newArray();
var nopdc5_c4=newArray();
var nopdc5_c5=newArray();
var nopdc5_c6=newArray();
var nopdc5_c7=newArray();
var nbpdc10=newArray();
var pdc10_c1=newArray();
var pdc10_c2=newArray();
var pdc10_c3=newArray();
var pdc10_c4=newArray();
var pdc10_c5=newArray();
var pdc10_c6=newArray();
var pdc10_c7=newArray();
var nbnopdc10=newArray();
var nopdc10_c1=newArray();
var nopdc10_c2=newArray();
var nopdc10_c3=newArray();
var nopdc10_c4=newArray();
var nopdc10_c5=newArray();
var nopdc10_c6=newArray();
var nopdc10_c7=newArray();


//---INITIALIZATION
min_size=1580; //Minimum cell size in pixels = 170 um²
max_size=6060; //Minimum cell size in pixels = 650 um²
//Name of the channels and various combinations of 
c1="Neon Green"; c2="Red/dsRNA"; c3="FarRed/spike";
c4=c1+" & "+c2;c5=c1+" & "+c3;c6=c2+" & "+c3;
c7=c1+" & "+c2+" & "+c3;
c8="pdc(5um)";
c9="pdc(10um)";
fluoth=66; //Fluorescence threshold for cell detection
minth=2; //Minimum fluorescence per channel to be counted as positive
pdcth=2; //Threshold to detect pdc cells
nbc=5; //Number of channels, for final measurements
dir = getDirectory( "Choose the Directory" );
list = getFileList(dir);

setBatchMode(true);
for (i=0; i<list.length; i++){
	ext=substring(list[i],lengthOf(list[i])-4,lengthOf(list[i]));
	if (ext==".czi"){//Only opens .czi files
		run("Bio-Formats (Windowless)", "open=[" +dir+list[i]+"]");
		roiManager("Reset");
		if (nSlices==nbc)main(list[i],i); //Only analyzes 5-channel stacks
		close();
	}
} 

close("ROI Manager");
//Prints resulting cell numbers for each category, in 4 distinct files: pdc_5, nopdc_5, pdc10, nopdc10
print("Name [PDC 5 um]\t pdc \t #cells \t #"+c1+" only \t #"+c2+" only \t #"+c3 + " \t #"+c4+"\t #"+c5+" \t #"+c6 + " \t #"+c7);
for (i=0;i<tot.length;i++){
	if (tot[i]>0){
		print(list[i]+ "\t w/ pdc \t" +nbpdc5[i]+" \t" +pdc5_c1[i]+" \t" +pdc5_c2[i]+" \t" +pdc5_c3[i]+" \t" +pdc5_c4[i]+" \t" +pdc5_c5[i]+" \t" +pdc5_c6[i]+" \t" +pdc5_c7[i]);
	}
}
selectWindow("Log");
saveAs("Text",dir+"results_pdc_5.txt");
run("Close");
print("Name [no PDC 5 um]\t pdc \t #cells \t #"+c1+" only \t #"+c2+" only \t #"+c3 + " \t #"+c4+"\t #"+c5+" \t #"+c6 + " \t #"+c7);
for (i=0;i<tot.length;i++){
	if (tot[i]>0){
		print(list[i]+ "\t no pdc \t" +nbnopdc5[i]+" \t" +nopdc5_c1[i]+" \t" +nopdc5_c2[i]+" \t" +nopdc5_c3[i]+" \t" +nopdc5_c4[i]+" \t" +nopdc5_c5[i]+" \t" +nopdc5_c6[i]+" \t" +nopdc5_c7[i]);
	}
}
selectWindow("Log");
saveAs("Text",dir+"results_nopdc_5.txt");
run("Close");

print("Name [PDC 10 um]\t pdc \t #cells \t #"+c1+" only \t #"+c2+" only \t #"+c3 + " \t #"+c4+"\t #"+c5+" \t #"+c6 + " \t #"+c7);
for (i=0;i<tot.length;i++){
	if (tot[i]>0){
		print(list[i]+ "\t w/ pdc \t" +nbpdc10[i]+" \t" +pdc10_c1[i]+" \t" +pdc10_c2[i]+" \t" +pdc10_c3[i]+" \t" +pdc10_c4[i]+" \t" +pdc10_c5[i]+" \t" +pdc10_c6[i]+" \t" +pdc10_c7[i]);
	}
}
selectWindow("Log");
saveAs("Text",dir+"results_pdc_10.txt");
run("Close");
print("Name [no PDC 5 um]\t pdc \t #cells \t #"+c1+" only \t #"+c2+" only \t #"+c3 + " \t #"+c4+"\t #"+c5+" \t #"+c6 + " \t #"+c7);
for (i=0;i<tot.length;i++){
	if (tot[i]>0){
		print(list[i]+ "\t no pdc \t" +nbnopdc10[i]+" \t" +nopdc10_c1[i]+" \t" +nopdc10_c2[i]+" \t" +nopdc10_c3[i]+" \t" +nopdc10_c4[i]+" \t" +nopdc10_c5[i]+" \t" +nopdc10_c6[i]+" \t" +nopdc10_c7[i]);
	}
}
selectWindow("Log");
saveAs("Text",dir+"results_nopdc_10.txt");
run("Close");
setBatchMode("exit and display");


function main(t,j){ //For each file named t and of index j...
	selectWindow(t);
	run("Z Project...", "start=1 stop=3 projection=[Max Intensity]");
	rename("temp");
	//... threshold the average projection of the three first channels...
	run("Smooth");
	setThreshold(fluoth, 65535, "raw");
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Fill Holes");
	run("Options...", "iterations=1 count=1 black do=Dilate");
	run("Options...", "iterations=2 count=1 black pad do=Erode");
	run("Options...", "iterations=1 count=1 black do=Dilate");
	run("Watershed");
	//... and keep the appropriate particles as cells
	run("Analyze Particles...", "size="+min_size+"-"+max_size+" pixel exclude show=Masks exclude");
	run("Grays");
	rename("cells"); close("temp");

	//XXX ATTENTION A LA VALIDATION COMPLEXE DES BORDS !!!
	//Builds a voronoi map from cell detection...
	run("Duplicate...", "title=voronoi");
	run("Voronoi");
	setAutoThreshold("Default dark");
	setThreshold(1, 255);
	setOption("BlackBackground", true);
	run("Convert to Mask");
	run("Invert");
	run("Analyze Particles...", "size="+min_size+"-Infinity pixel clear add");
	//... and keeps voronoi tiles in ROI Manager to retain same numbering for pdc detection (see below)
	nROIs=roiManager("Count");
	close("voronoi");

	//Builds a new "mask" stack that contains the original signal within the cells only (c1->c3) 
	//or within dilated cells (c4->c5) for pdc detection
	selectWindow("cells");
	run("Divide...", "value=255");
	imageCalculator("Multiply create stack", t,"cells");
	rename("mask");
	selectWindow("cells");
	run("Multiply...", "value=255.000");
	run("Duplicate...", "title=fat5");
	run("Options...", "iterations=15 count=1 black do=Dilate");
	run("Duplicate...", "title=fat10");
	run("Options...", "iterations=16 count=1 black do=Dilate");
	run("Divide...", "value=255");
	selectWindow("fat5");
	run("Divide...", "value=255");
	selectWindow(t);
	setSlice(4);
	imageCalculator("Multiply create", t,"fat5");
	run("Select All");
	run("Copy");
	close(); close("fat5");
	selectWindow("mask");
	setSlice(4);
	run("Paste");
	selectWindow(t);
	setSlice(4);
	imageCalculator("Multiply create", t,"fat10");
	run("Select All");
	run("Copy");
	close(); close("fat10");
	selectWindow("mask");
	setSlice(5);
	run("Paste");
	
	//Measures the detected cell areas within each voronoi tile
	selectWindow("cells");
	setAutoThreshold("Default dark");
	run("Set Measurements...", "area limit redirect=None decimal=3");
	roiManager("Measure");
	area=Table.getColumn("Area");
	close("Results");
	
	//Measures the mean fluorescence signal for each channel of the "mask" stack
	//within each voronoi tile
	means=newArray();
	means[0]=0;
	selectWindow("mask");
	run("Set Measurements...", "integrated redirect=None decimal=3");
	roiManager("Multi Measure");
	for (r=0;r<nROIs;r++){
		intdens=Table.getColumn("IntDen"+(r+1));
		means=Array.concat(means,intdens);
		for (c=0;c<3;c++){
			means[c+1+r*nbc]=means[c+1+r*nbc]/area[r];
		}
		//Estimate the area of the dilated cells to calculate the mean of channels 4 and 5
		means[5+r*nbc]=means[5+r*nbc]-means[4+r*nbc];
		R=sqrt(area[r]/PI)+5; //Estimated radius + 5 um for pdc5 = channel4
		area_fat5=PI*R*R;
		means[4+r*nbc]=means[4+r*nbc]/area_fat5;
		R=R+5;//Previously estimated radius + 5 um for pdc10 = channel5
		area_fat10=PI*R*R-area_fat5; //Considers only the annulus around previously calculated area_fat5
		means[5+r*nbc]=means[5+r*nbc]/area_fat10;
	}
	close("Results");

	//Calculates the distribution of variously labelled cells, for pdc5, nopdc5, pdc10 and nopdc10 cells
	tot[j]=nROIs; nbpdc5[j]=0; nbnopdc5[j]=0; nbpdc10[j]=0; nbnopdc10[j]=0;
	pdc5_c1[j]=0; pdc5_c2[j]=0; pdc5_c3[j]=0; pdc5_c4[j]=0; pdc5_c5[j]=0; pdc5_c6[j]=0; pdc5_c7[j]=0;
	nopdc5_c1[j]=0; nopdc5_c2[j]=0; nopdc5_c3[j]=0; nopdc5_c4[j]=0; nopdc5_c5[j]=0; nopdc5_c6[j]=0; nopdc5_c7[j]=0;
	pdc10_c1[j]=0; pdc10_c2[j]=0; pdc10_c3[j]=0; pdc10_c4[j]=0; pdc10_c5[j]=0; pdc10_c6[j]=0; pdc10_c7[j]=0;
	nopdc10_c1[j]=0; nopdc10_c2[j]=0; nopdc10_c3[j]=0; nopdc10_c4[j]=0; nopdc10_c5[j]=0; nopdc10_c6[j]=0; nopdc10_c7[j]=0;
	print("Name \t Area \t "+c1+" \t "+c2+" \t "+c3+" \t "+c8+" \t "+c9);
	for (r=0;r<nROIs;r++){
		roiManager("Select",r);
		n="cell-"+(r+1);
		roiManager("Rename", n);
		//Prints area and intensity profile of each cell
		print(n+" \t"+area[r]+ "\t"+means[1+r*nbc]+ "\t"+means[2+r*nbc]+ "\t"+means[3+r*nbc]+ "\t"+means[4+r*nbc]+ "\t"+means[5+r*nbc]);
		c1pos=false; c2pos=false; c3pos=false; c8pos=false; c9pos=false;
		if (means[1+r*nbc]>minth) c1pos=true; 
		if (means[2+r*nbc]>minth) c2pos=true;
		if (means[3+r*nbc]>minth) c3pos=true;
		if (means[4+r*nbc]>pdcth) c8pos=true;
		if ((means[4+r*nbc]>pdcth)|| (means[5+r*nbc]>pdcth)) c9pos=true;
		
		if (c8pos) { //c8pos detects whether there are pdc within 5 um of the scrutinized cell
			nbpdc5[j]++;
			if (c1pos & !c2pos & !c3pos) pdc5_c1[j]++;
			if (!c1pos & c2pos & !c3pos) pdc5_c2[j]++;
			if (!c1pos & !c2pos & c3pos) pdc5_c3[j]++;
			if (c1pos & c2pos & !c3pos) pdc5_c4[j]++;
			if (c1pos & !c2pos & c3pos) pdc5_c5[j]++;
			if (!c1pos & c2pos & c3pos) pdc5_c6[j]++;
			if (c1pos & c2pos & c3pos) pdc5_c7[j]++;
			}else{
			nbnopdc5[j]++;
			if (c1pos & !c2pos & !c3pos) nopdc5_c1[j]++;
			if (!c1pos & c2pos & !c3pos) nopdc5_c2[j]++;
			if (!c1pos & !c2pos & c3pos) nopdc5_c3[j]++;
			if (c1pos & c2pos & !c3pos) nopdc5_c4[j]++;
			if (c1pos & !c2pos & c3pos) nopdc5_c5[j]++;
			if (!c1pos & c2pos & c3pos) nopdc5_c6[j]++;
			if (c1pos & c2pos & c3pos) nopdc5_c7[j]++;
		}

		if (c9pos) { //c9pos detects whether there are pdc within 10 um of the scrutinized cell
			nbpdc10[j]++;
			if (c1pos & !c2pos & !c3pos) pdc10_c1[j]++;
			if (!c1pos & c2pos & !c3pos) pdc10_c2[j]++;
			if (!c1pos & !c2pos & c3pos) pdc10_c3[j]++;
			if (c1pos & c2pos & !c3pos) pdc10_c4[j]++;
			if (c1pos & !c2pos & c3pos) pdc10_c5[j]++;
			if (!c1pos & c2pos & c3pos) pdc10_c6[j]++;
			if (c1pos & c2pos & c3pos) pdc10_c7[j]++;
			}else{
			nbnopdc10[j]++;
			if (c1pos & !c2pos & !c3pos) nopdc10_c1[j]++;
			if (!c1pos & c2pos & !c3pos) nopdc10_c2[j]++;
			if (!c1pos & !c2pos & c3pos) nopdc10_c3[j]++;
			if (c1pos & c2pos & !c3pos) nopdc10_c4[j]++;
			if (c1pos & !c2pos & c3pos) nopdc10_c5[j]++;
			if (!c1pos & c2pos & c3pos) nopdc10_c6[j]++;
			if (c1pos & c2pos & c3pos) nopdc10_c7[j]++;
		}
	close("mask");
	}

	//Transforms voronoi tiles in cell contours before saving them
	selectWindow("cells");
	for (r=0;r<nROIs;r++){
		roiManager("Select", 0);
		n=Roi.getName();
		run("Analyze Particles...", "size="+min_size+"-Infinity pixel add");
		roiManager("Select", nROIs);
		roiManager("Rename", n);
		//Color the ROIs of cells closer than 5 um to PDC [~c8pos) in cyan
		//and those closer than 10 um to PDC (~c9pos) in magenta
		if (means[4+r*nbc]>pdcth){
			Roi.setStrokeColor("cyan");
			}else{if (means[5+r*nbc]>pdcth) Roi.setStrokeColor("magenta");
		}
		roiManager("Select", 0);
		roiManager("Delete");
	}
	close("cells"); 
	name=substring(t,0,lengthOf(t)-4);
	roiManager("Save",dir+name+".zip");
	selectWindow("Log");
	saveAs("Text",dir+name+".txt");
	close("Log");
}


