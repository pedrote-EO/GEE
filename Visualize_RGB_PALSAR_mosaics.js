/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Written by: P. Rodriguez-Veiga
//Purpose: Visualize PALSAR/PALSAR-2 mosaics 
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1 - INITIAL SETTINGS (Input needed)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set initial and final year for the annual mosaics to generate RGB. Set the same year for both if you only want 1 year output, or different for temporal composite
// PALSAR/PALSAR-2 mosaics available: 2007,2008,2009,2010,2015,2016,2017
var YiniDate1 = 2007;
var YendDate1 = 2008;

var YiniDate2 = 2010;
var YendDate2 = 2015;

var YiniDate3 = 2016;
var YendDate3 = 2017;

//Select date for single date RGB visualization
var dateRGB = '3'  //Date '1', '2' or '3' (as above)

//Select reducer type (mean, median, etc)
var s_reducer = ee.Reducer.mean()

Map.setOptions('SATELLITE');
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2 - FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for ALOS collection qa, and data transformations
var maskQA = function(image) {
   return image.mask(image.select('qa').neq(0)
   .and(image.select('qa').neq(100)
   .and(image.select('qa').neq(150)
   )));
  };

var to_dB = function(image) {
  // Keep this list of properties.
  var keepProperties = ['system:asset_size', 'system:footprint', 'system:index','system:time_start'];
  // apply function
  var log = image.multiply(image).log10().multiply(10).subtract(83);
  // Return a new Feature, copying properties from the old Feature.
  return ee.Image(log).copyProperties(image, keepProperties);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3 - APPLYING FUNCTIONS TO COLLECTION
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Filter temporary collection
var PALSARtemp = ee.ImageCollection('JAXA/ALOS/PALSAR/YEARLY/SAR')
        .map(maskQA)
        .map(to_dB);

//Select only polarization bands
var PALSAR = PALSARtemp.select(['HH','HV']);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 4 - TEMPORAL COMPOSITING
/////////////////////////////////////////////////////////////////////////////////////////////////////////
var iniDate1 = YiniDate1+'-01-01';
var endDate1 = YendDate1+'-12-31';
var iniDate2 = YiniDate2+'-01-01';
var endDate2 = YendDate2+'-12-31';
var iniDate3 = YiniDate3+'-01-01';
var endDate3 = YendDate3+'-12-31';

///Composites selected dates
var P_data1 = PALSAR
.filterDate(iniDate1, endDate1)
.reduce(s_reducer)
var P_data1 = P_data1.addBands(P_data1.select(0).divide(P_data1.select(1))).rename('HH1','HV1','ratio1');

var P_data2 = PALSAR
.filterDate(iniDate2, endDate2)
.reduce(s_reducer);
var P_data2 = P_data2.addBands(P_data2.select(0).divide(P_data2.select(1))).rename('HH2','HV2','ratio2');

var P_data3 = PALSAR
.filterDate(iniDate3, endDate3)
.reduce(s_reducer);
var P_data3 = P_data3.addBands(P_data3.select(0).divide(P_data3.select(1))).rename('HH3','HV3','ratio3');

/// Add all scenes together
var newALOS = P_data1.addBands(P_data2).addBands(P_data3);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 5 - VISUALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vis parameters
var imageVisParam = {"opacity":1,"bands":["HV1","HV2","HV3"],"max":-10, "min":-25,"gamma":0.75};
var imageVisParam1 = {"opacity":1,"bands":["HH"+dateRGB,"HV"+dateRGB,"ratio"+dateRGB],"max":[0,-10,1], "min":[-25,-25,0],"gamma":0.75};

Map.addLayer(newALOS, imageVisParam,'PALSAR_timeRGB');
Map.addLayer(newALOS, imageVisParam1,'PALSAR_RGB');
