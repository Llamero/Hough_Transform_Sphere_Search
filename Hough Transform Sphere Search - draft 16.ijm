//Set measurement tool to measure XY location and min/max of all objects in slice (Shape descriptors will be used later)
run("Set Measurements...", "area centroid perimeter fit shape feret's redirect=None decimal=9");
close("*");
if (isOpen("Results")) { 
    selectWindow("Results"); 
    run("Close"); 
}  
print("\\Clear");

setBatchMode(true);

//Sets the default lif file structure
dataSorting = "By collection within *.lif file";
//Sets the default as to record the hough scores
recordHough = true;
//sets the default to delete all intermediate files
deleteFiles = true;

//How much median filter to apply to the raw data
medianFilter = 3;
//How much to blur the sphere connection mask to allow fo a smooth expansion of the connections to be subtracted
connectionBlur = 1;
//How much to expand the saturated redion of connections (1 = total expansion, 255 = no expansion); - originally 125
connectionMax = 125;
//What Fraction of the stack (1/n) should be searched for the slide slice
slideFraction = 2;
//What factor should the mean stack intensity should be divided by to set the intensity threshold to start finding the slide.
findSlideFactor = 3;

//The lower threshold for a Hough score
houghThreshold = 40;

//Search radius limit when finding the median radius
minRadiusFactor = 0.8;
maxRadiusFactor = 1.2;

//Minimum and maximum ratio of median radius/Hough radius
minRadiusRatio = 0.88;
maxRadiusRatio = 1.15;

//The number of test lines used to characterized the x,y radii of ellipsoids in the cap search - must be a multiple of 2!
testLineCounter = 64;
//The percetile radius found to be used as the actual XY radius (0-1);
percentileRadius = 0.5;
//The ratio of the area of the z-projection of the image contained in the sphere to the cross-sectional area of the estimated sphere
areaRatioThreshold = 0.9;

//Minimum Hough search radius
houghMin = 10;
//Maximum Hough Search radius
houghMax = 150;

//Number of Hough circles to find per search
nCircles = 10;

//Ask user to choose the input and output directories
directory = getDirectory("Choose input directory");
rawFileList = getFileList(directory);
outputDirectory = getDirectory("Choose output directory");
run("Bio-Formats Macro Extensions");

useDefault = getBoolean("Do you want to use the default settings?");

if(!useDefault){
	//Create a user menu to over-ride default settings
	Dialog.create("Please enter values for the search algorithm:"); 
	Dialog.addNumber("3D median filter used to process raw stack:", medianFilter);
	Dialog.addNumber("Blur connections to smooth mask before subtraction:", connectionBlur);
	Dialog.addNumber("Threshold used to expand connections (1 = total expansion, 255 = no expansion):;", connectionMax);
	Dialog.addNumber("Fraction of the stack (1/n) that will be searched for the slide slice:", slideFraction);
	Dialog.addNumber("Fraction of the mean stack intensity (1/n) used to start finding the slide slice:", findSlideFactor);
	Dialog.addNumber("The minimum Hough score cutoff:", houghThreshold);
	Dialog.addNumber("Minimum inner limit used when searching for median radius (multiple of Hough radius):", minRadiusFactor);
	Dialog.addNumber("Maximum outer limit used when searching for median radius (multiple of Hough radius):", maxRadiusFactor);
	Dialog.addNumber("Minimum ratio of median radius / Hough radius:", minRadiusRatio);
	Dialog.addNumber("Maximum ratio of median radius / Hough radius::", maxRadiusRatio);
	Dialog.addNumber("Number of radial test lines used to find the median radius:", testLineCounter);
	Dialog.addNumber("Percentile of measured radii to be used as the actual radius (0.5 = median):", percentileRadius);
	Dialog.addNumber("Cutoff ratio of area of z-projection of mask : area of estimated sphere (higher = more stringent):", areaRatioThreshold);
	Dialog.addNumber("Minimum pixel radius to be used in the Hough transform search:", houghMin);
	Dialog.addNumber("Maximum pixel radius to be used in the Hough transform search:", houghMax);
	Dialog.addNumber("Number of circles to be returned per Hough transform iteration:", nCircles);
	Dialog.addCheckbox("Record the Hough transform score (rarely may cause search to fail):", true);
	Dialog.addCheckbox("Delete all intermediate files used in image processing:", true);	
	items = newArray("By *.lif file", "By collection within *.lif file");
  	Dialog.addRadioButtonGroup("How is the sample data organized:", items, 1, 3, "By collection within *.lif file");
	Dialog.show();
	
	//Save state of checked boxes
	medianFilter = Dialog.getNumber();
	connectionBlur = Dialog.getNumber();
	connectionMax = Dialog.getNumber();
	slideFraction = Dialog.getNumber();
	findSlideFactor = Dialog.getNumber();
	houghThreshold = Dialog.getNumber();
	minRadiusFactor = Dialog.getNumber();
	maxRadiusFactor = Dialog.getNumber();
	minRadiusRatio = Dialog.getNumber();
	maxRadiusRatio = Dialog.getNumber();
	testLineCounter = Dialog.getNumber();
	percentileRadius = Dialog.getNumber();
	areaRatioThreshold = Dialog.getNumber();
	houghMin = Dialog.getNumber();
	houghMax = Dialog.getNumber();
	nCircles = Dialog.getNumber();
	recordHough = Dialog.getCheckbox();
	deleteFiles = Dialog.getCheckbox();
	dataSorting = Dialog.getRadioButton();
	  
}

//Count how many files in the directory end in *.lif
lifCount = 0;
for(a=0; a<rawFileList.length; a++){
	if(endsWith(rawFileList[a], ".lif")){
		lifCount = lifCount + 1;
	}
}

//Build a new array of just the *.lif files in the directory
fileList = newArray(lifCount);
lifCount = 0;
for(a=0; a<rawFileList.length; a++){
	if(endsWith(rawFileList[a], ".lif")){
		fileList[lifCount] = rawFileList[a];
		lifCount = lifCount + 1;
	}
}
print("Processing the following files:");
Array.print(fileList);

//Count the number of images there are to be processed within the directory
nSamples = 0;
for (i=0; i<fileList.length; i++) {
	file = directory + fileList[i];
	Ext.setId(file);
	
	//Measure number of series in file
	Ext.getSeriesCount(nStacks);
	nSamples = nSamples + nStacks;
}

//Create an array for saving all sample sames within the directory
sampleNames = newArray(nSamples);
	
//Build an array of all sample names within the directory
sampleCounter = 0;
excludeCounter = 0;
for (i=0; i<fileList.length; i++) {
	file = directory + fileList[i];
	Ext.setId(file);
	
	//Measure number of series in file
	Ext.getSeriesCount(nStacks);

	//Open all stacks from set of lif files, one stack at a time
	for(j=0; j<nStacks; j++) {	

		//Only open the image if it has one / (i.e. part of collection, but not subcollection), and ends with a series number (not a, b, merged, etc.)
		Ext.setSeries(j);
		Ext.getSeriesName(seriesTitle);

		if(dataSorting == "By *.lif file"){
			//Make sure the series is not part of a collection
			if(!matches(seriesTitle, ".*/.*$") && matches(seriesTitle, ".*Series[0-9][0-9][0-9]$")){
				sampleNames[sampleCounter] = replace(fileList[i], "\\.lif.*", "");
			}
			//If this series name is not formatted correctly, exclude it from the list
			else{
				sampleNames[sampleCounter] = "N/A";
				excludeCounter = excludeCounter + 1;
			}
		}
		else{
			//Make sure the series is part of a collection, but not a sub-collection
			if(matches(seriesTitle, ".*/.*$") && !matches(seriesTitle, ".*/.*/.*$") && matches(seriesTitle, ".*Series[0-9][0-9][0-9]$")){
				seriesTitle = replace(seriesTitle, fileList[i] + " - ", "");
				sampleNames[sampleCounter] = replace(seriesTitle, "/.*", "");
			}
			//If this series name is not formatted correctly, exclude it from the list
			else{
				sampleNames[sampleCounter] = "N/A";
				excludeCounter = excludeCounter + 1;
			}	
		}

		//Count up one sample index
		sampleCounter = sampleCounter + 1;
	}
}

print(excludeCounter + " images/stacks were excluded from analysis due to incorrect naming format.");

//Record the number of series within the same sample:
sampleNamesCopy = Array.copy(sampleNames);
seriesNumbers = newArray(sampleNames.length);
Array.fill(seriesNumbers, 0);
//Search for a unique sample name
for(a=0; a<sampleNames.length; a++){
	seriesIDcounter = 1;
	
	//If found, replace name with "N/A" and record series number as 1
	if(sampleNamesCopy[a] != "N/A"){
		currentName = sampleNamesCopy[a];
		sampleNamesCopy[a] = "N/A";
		seriesNumbers[a] = seriesIDcounter;
		seriesIDcounter = seriesIDcounter + 1;

		//Search for all identical sample names and replace with N/A and record series number accordingly
		for(b=a; b<sampleNames.length; b++){
			if(sampleNamesCopy[b] == currentName){
				sampleNamesCopy[b] = "N/A";
				seriesNumbers[b] = seriesIDcounter;
				seriesIDcounter = seriesIDcounter + 1;
			}
		}
	}
}

//Create an array for tracking circle counts - index0 = # Found, index1 = # Overlap, index2 = # bad, index3 = boolean - good sphere found this iteration, index4 = iteration counter,
//index5 = # below Hough threshold, index6 = boolean - start removing bad spheres, index7 = boolean - good sphere found this cycle, index8 = total # found, index9 - # outside radius ratio 
circleTracker = newArray(10);

//Create an array for saving the names of all intermediate file names
processedImages = newArray(5);

//Open all data contained within the directory
sampleCounter = 0;

for (i=0; i<fileList.length; i++) {
	file = directory + fileList[i];
	Ext.setId(file);

	//Measure number of series in file
	Ext.getSeriesCount(nStacks);
	
	//Open all stacks from set of lif files, one stack at a time
	for(j=0; j<nStacks; j++) {	
		
		//Only open the image if it was formatted correctly (i.e. not listed as "N/A")
		if(sampleNames[sampleCounter] != "N/A"){
			run("Bio-Formats Importer", "open=file color_mode=Default split_channels view=[Standard ImageJ] stack_order=Default series_"+d2s(j+1,0));

			//Close the brightfield channel
			close("*C=1");
			
			//Retrieve stack absolute scaling information
			getVoxelSize(voxelWidth, voxelHeight, voxelDepth, voxelUnit);
			imageTitle = getTitle();

			//Create the final series name based off of the naming arrays
			seriesName = sampleNames[sampleCounter] + " - Series " + seriesNumbers[sampleCounter];
			print("Processing image " + sampleCounter  + 1 + " of " + nSamples + ": " + seriesName);

			//Reset the processed images array
			Array.fill(processedImages,0);
			retryCounter = 0;

			//Find the slide slice and then create a mask		
			slideSlice = findSlide(imageTitle);
			processedImages = makeBinaryMask(imageTitle, processedImages);
			binaryMask = processedImages[4];
			
			//Create a new image stack to store the found labels
			selectWindow(binaryMask);
			getDimensions(stackWidth, stackHeight, dummy, stackSlices, dummy);
			labels = seriesName + " - sphere labels.tif";
			newImage(labels, "8-bit black", stackWidth, stackHeight, stackSlices);
			
			//Reset the circle tracker array
			Array.fill(circleTracker, 0);

			//Search for all spheres in the mask
			print("Starting sphere search...");
			while(circleTracker[5]<nCircles || !circleTracker[6]){
			
				//Count up one on the iteration counter
				circleTracker[0] = 0;
				circleTracker[1] = 0;
				circleTracker[2] = 0;
				circleTracker[3] = 0;
				circleTracker[4] = circleTracker[4] + 1;
				circleTracker[9] = 0;

				//Use the hough threshold counter to decide to stop search, if threshold is available
				if(recordHough && houghThreshold >= 1){
					circleTracker[5] = 0;
				}
				//Otherwise, exclude hough threshold from search criteria
				else{
					circleTracker[5] = nCircles;
				}
			
				houghSearch(binaryMask);
				for(circle = 1; circle<=nCircles; circle++){
					//Reset whether a good sphere was found to false
					circleTracker[7] = 0;
					circleTracker = findHoughSphere(binaryMask, labels, circleTracker);
				}
				close("Centroid Map contains " + nCircles + " circles");

				//Report Current Progress to the log
				print("Iteration #" + circleTracker[4] + " Round #" + circleTracker[6] + 1 + ": " + circleTracker[0] + " approved, " + circleTracker[1] + " rejected - overlapping centroid, " + circleTracker[2] + " rejected - poor quality, " + circleTracker[5] + " rejected - low Hough score, " + circleTracker[9] + " rejected - outside radius range.");

				//If the "remove only good spheres" search is exhausted, then begin removing bad spheres too
				if(!circleTracker[3] && !circleTracker[6]){
					circleTracker[6] = 1;
				}

				//If it is round 2 and the retry failed to produce a new sphere, count up one on the retry counter
				if(!circleTracker[3] && circleTracker[6]){
					retryCounter = retryCounter + 1;
				}
				//Otherwise, reset the retry counter
				else{
					retryCounter = 0;
				}

				//If circle search is complete, then save results and close the results window
				if(!circleTracker[3] && (circleTracker[5] == nCircles || retryCounter > 3) && circleTracker[6]){
					//Set circleTracker[5] to nCircles to ensure the while loop is exited in case the search stops due to retryCounter instead
					circleTracker[5] = nCircles;
					print("Search complete!");

					//Save the results table
					saveAs("Results", outputDirectory + seriesName + " - measurements.xls");

					//Open the 3x3x3 median file and merge it with the labels file and save the composite
					open(outputDirectory + processedImages[0]);
					run("Merge Channels...", "c1=[" + labels + "] c2=[" + processedImages[0] + "] create ignore");
					saveAs("Tiff", outputDirectory + seriesName + " - median and labels composite.tif");
					close("*");
 				    selectWindow("Results"); 
				    run("Close");

					//If selected, delete all intermediate images
					//don't delete the first image in the array, as it is the 3x3x3 median and is worth keeping
					if(deleteFiles){
						for(a=0; a<processedImages.length; a++){
							dummy = File.delete(outputDirectory + processedImages[a]);
						}
					}
				}
			}
		}
		//Increment up the sample counter
		sampleCounter = sampleCounter + 1;	
	}
}

setBatchMode(false);

//Save the processing log
selectWindow("Log");
getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
logName = "" + year + "-" + month + "-" + dayOfMonth + " sphere search log at " + hour + " hours " + minute + " minutes";
saveAs("Text", outputDirectory + logName);

function findSlide(imageTitle){
	showStatus("Finding the slide slice...");
	selectWindow(imageTitle);
	
	//Find the brightest slice in the stack, as this will be the slide (to be used as reference in analysis)
	//Initialize the mean intensity measurement and slide slice # variable
	slideSlice = 0;
	Stack.getStatistics(dummy, stackMean, dummy, stackMax, dummy);

	//Start the maxMean (slide) at the max intensity in the whole stack, so that search doesn't start searching for the slide at the very start of the search
	maxMean = stackMax;
	startSearch = stackMean/findSlideFactor; 
	
	for (b=1; b<=nSlices/slideFraction; b++){
		setSlice(b);
		getStatistics(dummy, mean);

		//if you find a slice that is brighter than the starting cutoff, start searching for the slide by setting max to 0
		if (mean>startSearch && slideSlice == 0){
			maxMean = 0;			
		}
		if (mean>maxMean){
			slideSlice = b;
			maxMean = mean;
		}
		//If you are past the first maximmum, stop searching for the slide
		if (mean < maxMean && slideSlice > 0){
			b = nSlices + 1;
		}
	}
	return slideSlice;
}

//Create a binary mask of the raw data
function makeBinaryMask(imageTitle, processedImages){
	selectWindow(imageTitle);
	//------------------------------------Remove the bright connection points between the spheres--------------------------------------------------------------------------------------
	//Where two spheres touch there is a disproportionately bright spot.  When using intensity in Z as a reference to find spheres, these connections can be towards the top of the spheres, and
	//therefore will result in the false identification of spheres (i.e. a sphere will be found centered on the connection).  Therefore, the connections need to be removed first, to generate an image of
	//just the spheres.

	//Remove noise with 3D median
	run("8-bit");
	run("Median 3D...", "x=" + medianFilter + " y=" + medianFilter + " z=" + medianFilter);

	//Save the filtered data
	imageTitle = seriesName + " - " + medianFilter + "x" + medianFilter + "x" + medianFilter + " median.tif";
	processedImages[0] = imageTitle;
	saveAs("Tiff", outputDirectory + imageTitle);
	
	//Autothreshold based on the maximum entropy of the stack histogram (proved very specific autothreshold for connections), and make a corresponding mask
	setAutoThreshold("MaxEntropy dark stack");
	run("Convert to Mask", "method=Default background=Default black");
	
	//Expand out the binarized connections by performing a Gaussian blur and then autocontrasting to re-saturate the center of each connection for complete removal
	run("Gaussian Blur 3D...", "x=" + connectionBlur + " y=" + connectionBlur + " z=" + connectionBlur + "");
	setMinAndMax(0, connectionMax);
	run("Apply LUT", "stack");

	
	//Save the connections image to distinguish it form the original for the image calculator
	connectionTitle = replace(imageTitle, ".tif", " - sphere connections.tif");
	processedImages[1] = connectionTitle;
	saveAs("Tiff", outputDirectory + connectionTitle);
	
	//Open the 8-bit median filtered image and subtract the connections from it
	open(outputDirectory + imageTitle);
	imageCalculator("Subtract create stack",  imageTitle, connectionTitle);
	
	//Close the original image and the connections image
	close(imageTitle);
	close(connectionTitle);

	//Save the new image with the connections removed
	selectWindow("Result of " + imageTitle);
	ConnectionsCropped = replace(imageTitle, ".tif", " - connections cropped.tif");
	processedImages[2] = ConnectionsCropped;
	saveAs("Tiff", outputDirectory + ConnectionsCropped);

	//Autothreshold the stack and binarize
	selectWindow(ConnectionsCropped);
	setAutoThreshold("Mean dark stack");
	run("Convert to Mask", "method=Default background=Default black");

	//Save the filtered and binarized autofluroescence image
	binaryTitle = replace(ConnectionsCropped, ".tif", " - binarized.tif");
	processedImages[3] = binaryTitle;
	saveAs("Tiff", outputDirectory + binaryTitle);

	//Also save the image as "spheres cropped" to allow for tracking cropping in a separate image
	binaryCropped = replace(binaryTitle, ".tif", " - spheres cropped search.tif");
	processedImages[4] = binaryCropped;
	saveAs("Tiff", outputDirectory + binaryCropped);

	
	return processedImages;
}


function houghSearch(binaryMask){
	selectWindow(binaryMask);

	//Create a sum projection to reduce stack dimensionality
	run("Z Project...", "projection=[Sum Slices]");
	selectWindow("SUM_" + binaryMask);

	//Fina edges and binarize for hough transform
	run("Find Edges");
	setAutoThreshold("MaxEntropy dark");
//	setAutoThreshold("RenyiEntropy dark");
//	setAutoThreshold("Yen dark");
	setOption("BlackBackground", false);
	run("Convert to Mask");
	run("Grays");

	//Perform Hough transform and find best circle
	if(recordHough){
		run("Hough Circles Edit4", "minimum=" + houghMin + " maximum=" + houghMax + " increment=2 number=" + nCircles + " threshold=60 map reverse record");
	}
	else{
		run("Hough Circles Edit4", "minimum=" + houghMin + " maximum=" + houghMax + " increment=2 number=" + nCircles + " threshold=60 map reverse");
	}
	close("SUM_" + binaryMask);
}

function findHoughSphere(binaryMask, labels, circleTracker){
	//Find the centroid and record the radius
	houghCentroid = "Centroid Map contains " + nCircles + " circles";
	selectWindow(houghCentroid);
	getDimensions(width, height, channels, slices, frames);

	for (x=0; x<width; x++){
		for (y=0; y<height; y++){
			radius = getPixel(x, y);
			if(radius > 0){
				xCentroid = x;
				yCentroid = y;

				//Erase this centoid from the map
				setPixel(x,y,0);
				
				//stop search
				x = width;
				y = height;
				
			}
		}
	}

	//Retrieve the corresponding Hough score if available
	if(recordHough){
		houghScore = getPixel(xCentroid + 1, yCentroid);
		setPixel(xCentroid + 1, yCentroid, 0);
	}
	else{
		houghScore = 0;
	}

	

	//Find the corresponging z centroid in the binary mask
	//Start search at 1/4 measured radius above the slide slice (i.e. sphere cannot be sitting below the slide)
	selectWindow(binaryMask);
	startSearch = round(0.25 * radius + slideSlice);

	//Create an array to save all measured median radii
	medianRadii = newArray(nSlices+1);
	Array.fill(medianRadii, 0);

	//Create an array to save the measured radii
	sphereBoundary = newArray(testLineCounter);

	//Search through the entire stack for the largest median radius
	for(a=startSearch; a<nSlices; a++){
		setSlice(a);

		//Find the median radius of the sphere cross section at this slices
		//Prefill the array with min radii.  This prevents the array from returning zeros which crashes the search
		Array.fill(sphereBoundary, radius*maxRadiusFactor);
	
		//Starting at the deflated radius, measure out until you reach the boundary of the sphere, and record this radius
		for	(d=0; d<testLineCounter; d++){
			thetaRadians = (2 * PI) * (d/testLineCounter);
			//Start at the pixelRadius, because otherwise bubbles in the middles will cause false early terminations
			for(lineRadius=round(radius*minRadiusFactor); lineRadius<=(radius*maxRadiusFactor); lineRadius++){
				xLine = xCentroid + round(lineRadius*cos(thetaRadians));
				yLine = yCentroid + round(lineRadius*sin(thetaRadians));
				testPixel = getPixel(xLine, yLine);
		
				//If the test pixel is < 255, then record the radius, as this is the edge of the sphere in this direction
				if(testPixel < 255){
					sphereBoundary[d] = lineRadius;
					lineRadius = radius*maxRadiusFactor + 1;
				}
			}
		}
	
		//Find the desired percentile radius in the array.  Percentiles are used, as pits or fused spheres can cause large outliers
		Array.sort(sphereBoundary);
		medianRadii[a] = sphereBoundary[round(percentileRadius*testLineCounter)];			
	}

	//Convert the radius array to an image to find the peak radius
	newImage("median array", "8-bit black", medianRadii.length, 1, 1);
	selectWindow("median array");
	for(a=0; a<medianRadii.length; a++){
		setPixel(a,0,medianRadii[a]);
	}

	//Smooth the array so that local peaks are removed (remove high freq noise)
	selectWindow("median array");
	run("Gaussian Blur...", "sigma=3");

	//Find the maximum peak
	getStatistics(dummy, dummy, dummy, blurArrayMax);
	run("Find Maxima...", "noise=" + blurArrayMax + " output=[Point Selection]");
	getSelectionBounds(zCentroid, dummy, dummy, dummy);
	close("median array");

	//Check to make sure this centroid is not already occupied by another label
	selectWindow(labels);
	setSlice(zCentroid);
	labelOverlap = getPixel(xCentroid,yCentroid);
	
	//When creating the hypothetical sphere, use the Hough radius or measured median radius, whichever is greater
	if(radius>blurArrayMax){
		sphereRadius = radius;
	}
	else{
		sphereRadius = blurArrayMax;
	}
	
	//Create a new image stack to store the found labels
	selectWindow(binaryMask);
	getDimensions(stackWidth, stackHeight, dummy, stackSlices, dummy);
		
	//Draw the estimated sphere in a separate stack to allow for checking oringal object sphericity
	run("3D Draw Shape", "size=" + stackWidth + "," + stackHeight + "," + stackSlices + " center=" + xCentroid + "," + yCentroid + "," + zCentroid + " radius=" + sphereRadius + "," + sphereRadius + "," + sphereRadius + " vector1=1.0,0.0,0.0 vector2=0.0,1.0,0.0 res_xy=1 res_z=1 unit=pixel value=65535 display=[New stack]");

	selectWindow("Shape3D");
	run("8-bit");

//----------------------------------Crop the data from the original image contained within the approximate sphere, and check to make sure it too is spherical------------------------------------------
	//Speed up the divide by only dividing where the sphere is
	makeOval(xCentroid - (sphereRadius + 3), yCentroid - (sphereRadius + 3), 2*(sphereRadius + 3), 2*(sphereRadius + 3));
	run("Divide...", "value=255 stack");
	run("Select None");

	//If Hough score is above threshold, keep the sphere
	if(houghScore >= houghThreshold || houghScore == 0){

		//If the centroids do not overlap...
		if(labelOverlap == 0){

			//If the median radius is within range
			radiusRatio = blurArrayMax/radius;
			if(radiusRatio < maxRadiusRatio && radiusRatio > minRadiusRatio){

				//Multiply the original image by the estimated sphere to get the object in the original image predicted to be a sphere
				imageCalculator("Multiply create stack", binaryMask,"Shape3D");
				
				//Generate a mean projection of the result
				selectWindow("Result of " + binaryMask);
				run("Z Project...", "projection=[Average Intensity]");
				
				//Close the result form the image calculator as it is no longer needed
				close("Result of " + binaryMask);
				
				//Threshold the resulting projection for analysis
				//Originally Huang was used, but MaxEntropy allows better coverage of spheres that have large bubbles
				selectWindow("AVG_Result of " + binaryMask);					
				setAutoThreshold("Huang dark");
				run("Convert to Mask");
			
				//Measure the thresholded image
				//Sometimes the mask can result in small satellite particles which can show up in the results and therefore need
				//to be removed, otherwise, the satellite particle result may be chosen at random, which will not pass the 
				//quality filters and stop the sphere search prematurely.  Therefore, only particles above the min area ratio threshold
				//will be kept
				minimumParticleArea = (3.14159 * sphereRadius * sphereRadius * voxelWidth * voxelWidth)*areaRatioThreshold;
				selectWindow("AVG_Result of " + binaryMask);
				run("Analyze Particles...", "size=" + minimumParticleArea + "-Infinity circularity=0.00-1.00 display");
			
				//Close the average projection as it is no longer needed
				close("AVG_Result of " + binaryMask);
		
				//This needs to be conditional as otherwise if no measurement was made (area too small), then there would have been no result
				if(nResults > circleTracker[8]){
	
					//Record the roundness as a "quality score" and parameters for the given sphere
					setResult("X_centroid", circleTracker[8], xCentroid*voxelWidth);
					setResult("Y_centroid",circleTracker[8], yCentroid*voxelHeight);
					setResult("Z_centoid",circleTracker[8],zCentroid*voxelDepth);
					setResult("Hough_Radius",circleTracker[8],radius*voxelWidth);
					setResult("Median_Radius",circleTracker[8], blurArrayMax*voxelWidth);
					setResult("Hough_Score",circleTracker[8],houghScore);
					setResult("Method",circleTracker[8],"Hough");
					setResult("Unit_of_Measure",circleTracker[8],voxelUnit);
					updateResults();
				
					//Caclulate the ratio of the area of the actual sphere projection to the estimated sphere projection
					imageArea = getResult("Area",nResults-1);
					areaRatio = imageArea / (3.14159 * sphereRadius * sphereRadius * voxelWidth * voxelWidth);
					setResult("Area_Ratio",circleTracker[8],areaRatio);
					updateResults();
			
					//Increment the circle tracker found index one
					circleTracker[0] = circleTracker[0] + 1;
					circleTracker[8] = circleTracker[8] + 1;
			
					//Record that a circle has been found
					circleTracker[3] = 1;
					circleTracker[7] = 1;
			
					//Add the result to the labels stack
					//Adjust the sphere intensity to match the label count ID
					selectWindow("Shape3D");
					run("Multiply...", "value=" + circleTracker[8] + " stack");
					imageCalculator("Add stack", labels,"Shape3D");
				}
				//Otherwise, the sphere is not good quality so record that a poor quality sphere was rejected
				else{
					circleTracker[2] = circleTracker[2] + 1;
				}
			}
			else{
				//Increment the radius ratio tracker idex by one
				circleTracker[9] = circleTracker[9] + 1;
			}
		}
		else{
		//Increment the circle tracker overlap index one
		circleTracker[1] = circleTracker[1] + 1;
		}
	}
	else{

		//Increment the below Hough threshold index by one
		circleTracker[5] = circleTracker[5] + 1;
	}
	//Delete the sphere if it is only a good sphere on the first search pass, or remove all spheres (good or bad) on the second pass
	if((!circleTracker[6] && circleTracker[7]) || circleTracker[6]){

		//Speed up the multiple by only dividing where the sphere is
		selectWindow("Shape3D");
		makeOval(xCentroid - (sphereRadius + 3), yCentroid - (sphereRadius + 3), 2*(sphereRadius + 3), 2*(sphereRadius + 3));
		run("Multiply...", "value=255 stack");
		run("Select None");
		
		//Remove the found label from the mask and add it to the labels stack
		imageCalculator("Subtract stack", binaryMask,"Shape3D");
	}
	close("Shape3D");

	return circleTracker;
}


