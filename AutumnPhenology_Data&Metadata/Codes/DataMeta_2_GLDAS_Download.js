// This is the code for extract climate variables we need from GLDAS database on Google earth engine(GEE)


// load the GLDAS version 2.0 and filter the target years
var gldas_2_0 = ee.ImageCollection("NASA/GLDAS/V20/NOAH/G025/T3H").filterDate('1948-01-01','2001-01-01');
// load the GLDAS version 2.1 and filter the target years
var gldas_2_1 = ee.ImageCollection("NASA/GLDAS/V021/NOAH/G025/T3H").filterDate('2001-01-01','2019-01-01');
var gldas = gldas_2_0.merge(gldas_2_1);

// loadd the sample points information
var pointsToSample = ee.FeatureCollection("users/leonidmoore/AutaumPhenologyPoints");
Map.addLayer(pointsToSample, {}, "Points to Sample", false);
// print the inforamtion
// print(gldas)

var selectVariables = ['LWdown_f_tavg','SSWdown_f_tavg','Rainf_tavg','Tair_f_inst','Qair_f_inst']

for (var i = 0, len = selectVariables.length; i < len; i++) {
  // here words[i] is the array element
  // import the GLADIS data version 2 which ranges from 1948-2010;
  var gldas_2_Chose = gldas.select(selectVariables[i]);
  // print('Title 01:',gladis_2_Chose.first());

  var dataProjection = ee.Image(gldas.first()).projection();
  // print('Standard Projection',dataProjection);
  for (var year= 1948; year < 2019; year++){
    // Start and End Dates
    var inidate = ee.Date.fromYMD(year,1,1)
    var enddate = ee.Date.fromYMD((year+1),1,1)

    // Difference between start and end in days 
    var difdate = enddate.difference(inidate, 'day')
    // generate a list of the dates for everyday
    var lapse = ee.List.sequence(0, difdate.subtract(1))
    var inidate = ee.Date(year+'-01-01')
    var listDates = lapse.map(function(day){
      return inidate.advance(day, 'day')
    })
    // check the information of the listdates
    // print(listDates)

    // us map function to do this data reduction from 3- hourly to daily
    // Only the temperature data needs the maximum and minimum
    var composites = ee.ImageCollection.fromImages(listDates.map(function(d) {
      var filtered = gldas_2_Chose.filterDate(ee.Date(d).getRange('day'));
      // here is the place you have to change to .mean(), .max() and min()
      return filtered.mean().reproject({
        crs:dataProjection
      }).set('day', d);
    }));
    // .set('day', d).projection(composites)

    // print('Information of the first layer',composites.limit(1))
    // print('output',composites);
    // the code below is just to add the first layer "1948-01-01" to the google map in the UI.
    Map.addLayer(composites.limit(1), {}, "Composites First Layer", false);

    // Create an empty image to fill
    var emptyImage = ee.Image([]);

    // Iterate through the collection to make the new multiband image
    var multibandImage = ee.Image(composites.iterate(function(image, result) {
	    return ee.Image(result).addBands(image.reproject({
        crs:dataProjection
      }));
    }, emptyImage));
    // print(multibandImage)

    Map.addLayer(multibandImage, {}, "Multi Band", false)
    print('multi band image',multibandImage)
    // Using the sample regions function
    var pointsExtractedVariables = multibandImage.sampleRegions({
    	collection: pointsToSample
    });


    Export.table.toDrive({
      collection: pointsExtractedVariables,
      description: selectVariables[i]+"_"+year+"_Min_Extracted",
      folder:"20190301_Autaum_Phenology_Export",
      fileFormat:"CSV"
    })
  }
}

