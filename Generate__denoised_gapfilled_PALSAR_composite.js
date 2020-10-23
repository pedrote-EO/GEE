/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Main code body written by: Pedro Rodriguez-Veiga
//Purpose: Generate ALOS PALSAR / ALOS-2 PALSAR-2 mosaics denoised using multi-temporal filtered, and gap-filled
//Current version: 1.0
//Multi-temporal speckle filter -> https://www.researchgate.net/publication/3202307_Multitemporal_ERS_SAR_analysis_applied_to_forest_mapping
//Implementation based on code written by Gennadii Donchyts and Anouk Ville (Version: 1.0)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 1 - INITIAL SETTINGS (Input needed)
/////////////////////////////////////////////////////////////////////////////////////////////////////////

// Set the year for PALSAR/PALSAR-2 annual mosaic. Set the same year for both if you only one 1 year output, or different for temporal composite
// PALSAR/PALSAR-2 mosaics available so far in GEE: 2007,2008,2009,2010,2015,2016,2017,(2018, 2019, and ahead should be available in the future)
var YiniDate = 2015;
var YendDate = 2015;

//Type of reducer for the temporal compositing
//NOTE: If you use annual composites it does not matter which one you choose
var reducer_type = ee.Reducer.mean(); //mean, min or median,

//Moving window size for the multi-temporal filter (radius = 5 means 5x5 pixels)
//It has to be an odd number (e.g.3,5,7,etc)
var radius = 3;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 2 - FUNCTIONS
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Functions for ALOS collection QA to mask no data, radar layover and radar shadowing
var maskQA = function(image) {
   return image.mask(image.select('qa').neq(0)
   .and(image.select('qa').neq(100)
   .and(image.select('qa').neq(150)
   )));
  };

//Transform from DN to dB
var to_dB = function(image) {
  // Keep this list of properties.
  var keepProperties = ['system:asset_size', 'system:footprint', 'system:index','system:time_start'];
  // apply function
  var log = image.multiply(image).log10().multiply(10).subtract(83);
  // Return a new Feature, copying properties from the old Feature.
  return ee.Image(log).copyProperties(image, keepProperties);
};

//Transform from dB to power 
var dB_to_power = function(image) {
  // Keep this list of properties.
  var keepProperties = ['system:asset_size', 'system:footprint', 'system:index','system:time_start'];
  // apply function
  var power = image.expression('10**(0.1*image)',{'image': image})
  // Return a new Feature, copying properties from the old Feature.
  return ee.Image(power).copyProperties(image, keepProperties);
};

//Transform from power to dB
var power_to_dB = function(image) {
  // Keep this list of properties.
  var keepProperties = ['system:asset_size', 'system:footprint', 'system:index','system:time_start'];
  // apply function
  var db = image.abs().log10().multiply(10.0);
  // Return a new Feature, copying properties from the old Feature.
  return ee.Image(db).copyProperties(image, keepProperties);
};


/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 3 - MULTITEPORAL SPECKLE FILTER function
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//image is the original image, images is the temporal collection of images

function multitemporalDespeckle(images, radius, units, opt_timeWindow) {
  var timeWindow = opt_timeWindow || { before: -11, after: 11, units: 'year' }
  
  var bandNames = ee.Image(images.first()).bandNames()
  var bandNamesMean = bandNames.map(function(b) { return ee.String(b).cat('_mean') })
  var bandNamesRatio = bandNames.map(function(b) { return ee.String(b).cat('_ratio') })
  
  // compute space-average for all images
  var meanSpace = images.map(function(i) {
    var reducer = ee.Reducer.mean()
    var kernel = ee.Kernel.square(radius, units)
    
    var mean = i.reduceNeighborhood(reducer, kernel).rename(bandNamesMean)
    var ratio = i.divide(mean).rename(bandNamesRatio)

    return i.addBands(mean).addBands(ratio)
  })

  /***
   * computes a multi-temporal despeckle function for a single image
   */
  function multitemporalDespeckleSingle(image) {
    var t = image.date()
    var from = t.advance(ee.Number(timeWindow.before), timeWindow.units)
    var to = t.advance(ee.Number(timeWindow.after), timeWindow.units)
    var keepProperties = ['system:asset_size', 'system:footprint', 'system:index','system:time_start'];

    var meanSpace2 = ee.ImageCollection(meanSpace).select(bandNamesRatio).filterDate(from, to)
  
    var b = image.select(bandNamesMean)

    return b.multiply(meanSpace2.sum()).divide(meanSpace2.count()).rename(bandNames).copyProperties(image, keepProperties)
  }
  
  return meanSpace.map(multitemporalDespeckleSingle).select(bandNames)
}

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 4 - FILTERING COLLECTIONS AND PROCESSING
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//Filter temporary collection for multi-temporal filter
var PALSARtemp = ee.ImageCollection('JAXA/ALOS/PALSAR/YEARLY/SAR')
        .filterDate('2006-01-01', '2024-12-31')
        //Apply  QA function to convert
        .map(maskQA)
        //Transform DN to dB
        .map(to_dB);

//Select only polarization bands and transform from dB to power
var PALSAR = PALSARtemp.select(['HH','HV']).map(dB_to_power);
print(PALSAR);

// Apply multi-temporal filter (denoise images). It has to be done in power units
var units = 'pixels';
var PALSARDenoised = multitemporalDespeckle(PALSAR, radius, units);
print(PALSARDenoised);

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 5 - TEMPORAL COMPOSITING
/////////////////////////////////////////////////////////////////////////////////////////////////////////
var iniDate = YiniDate+'-01-01';
var endDate = YendDate+'-12-31';


/// Composites without multi-temporal filter for the selected dates and converting back to dB
var Original_data = PALSAR.map(power_to_dB)
.filterDate(iniDate, endDate)
.reduce(reducer_type).rename('HH','HV');
print(Original_data);

///Composites with multi-temporal filter for the selected dates and converting back to dB
var Denoised_data = PALSARDenoised.map(power_to_dB)
.filterDate(iniDate, endDate)
.reduce(reducer_type).rename('HH','HV');
print(Denoised_data);

/// Rename variable
var PALSARd = Denoised_data;

/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 6 - GAP FILLING
/////////////////////////////////////////////////////////////////////////////////////////////////////////
//Focal mean moving window to fill gaps
var P2texture7 = PALSARd.focal_mean({
  radius:3,
  kernelType: 'circle',
  iterations:7
});

//Sequential gap filling process
var PALSARdg = PALSARd.unmask(P2texture7);
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// 7 - VISUALIZATION
/////////////////////////////////////////////////////////////////////////////////////////////////////////
// Vis parameters

Map.setOptions('SATELLITE');
var imageVisParam = {"opacity":1,"bands":["HV","HH","HV"],"max":-5, "min":-20,"gamma":0.5};
Map.addLayer(Original_data, imageVisParam,'PALSAR2_original');
Map.addLayer(PALSARd, imageVisParam,'PALSAR2_denoised');
