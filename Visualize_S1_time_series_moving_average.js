/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Main code body written by: Pedro Rodriguez-Veiga
//Purpose: S1 moving average time series plot
//link: https://code.earthengine.google.com/9c76c5c257fcbb60dc4454a145eb8d91
//Current version: 1.0
/////////////////////////////////////////////////////////////////////////////////
/* Setting variables */
var iniDate = '2018-01-01';
var endDate = '2019-12-31';

//Manually create 1 small polygon over your area of interest (geometry)
Map.centerObject(geometry, 14);
/////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
/* Convert to Gamma 0 function */
function toGamma0(image) {
  var to_VV = image.select(0).subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0));
  var to_VH = image.select(1).subtract(image.select('angle').multiply(Math.PI/180.0).cos().log10().multiply(10.0));
  var gamma = to_VV.addBands(to_VH).rename('VV','VH')
  return gamma}

///////////////////////////////////////////////////////////////////
/* Processing code */
//collection for change-detection
var microwv = ee.ImageCollection("COPERNICUS/S1_GRD")
  .filterBounds(geometry)
  .filterDate(iniDate,endDate)
  .filterMetadata('instrumentMode', 'equals', 'IW')
  .filterMetadata('transmitterReceiverPolarisation', 'equals', ['VV','VH'])
  // .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'))
    ;
print(microwv);

var im1 = ee.Image(microwv.filterDate(iniDate,endDate).map(toGamma0).median());
var im1 = im1.addBands(im1.select(0).divide(im1.select(1))).rename('VV','VH','ratio');

// Vis parameters
var imageVisParam = {"opacity":1,"bands":["VH","VV","ratio"],"max":[-10,-5,1], "min":[-25,-20,-1],"gamma":0.75};

Map.addLayer(im1, imageVisParam,'S1_RGBmean',true)

var movAvg = function(collection, dt){
  // moving average filter 
  // collection: colection to smooth 
  // dt: smooth window. Time that will be added forward and back
  // dt = 1 will give a [t-1, t, t+1] window
  // uses a join on the date
  var joinFilter = ee.Filter.and(
    ee.Filter.maxDifference(
      dt*1*1000*60*60*24, 'system:time_start', null, 'system:time_start'
    )
  );
  
  var joinResult = ee.Join.saveAll('matches').apply(collection, collection, joinFilter);
  print('joinResult', joinResult);
  
  var averageMatchedImages = function (imageWithMatches) {
    var matchCollection = ee.ImageCollection.fromImages(imageWithMatches.get('matches'));
    var imageMean = matchCollection.reduce(ee.Reducer.mean());
    return ee.Image(imageWithMatches).addBands(imageMean);
  };

  var windowedAvgCollection = joinResult.map(averageMatchedImages);
  print('avgCollection', windowedAvgCollection);
  
  return windowedAvgCollection;
};

// Average every X days (e.g. 30 days)
var avgS1 = movAvg(microwv.select(['VV','VH']),30);
print(avgS1);

////////////////////////////////////////////////////////////////////////////////
/* Chart functions */

// Create an image time series chart.
var chart = ui.Chart.image.series({
  imageCollection: avgS1,
  region: geometry,
  reducer: ee.Reducer.mean(),
  scale: 10,
//  xProperty:'month'
});
chart.setOptions({
  title: 'Sentinel 1',
  vAxis: {
    title: 'Gamma 0 Backscatter (dB)'
  },
  legend: 'band',
  lineWidth: 1,
});

// Add the chart to the map.
chart.style().set({
  position: 'bottom-right',
  width: '500px',
  height: '300px'
});
Map.add(chart);
