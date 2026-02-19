var aoi = ee.FeatureCollection('projects/ee-yudong9607/assets/paper/Buvuma_down2018');
var all_points = ee.FeatureCollection('projects/ee-yudong9607/assets/paper/Test_thin43_2')
  .map(function(f) {
    return f.set('location_id', f.get('down_id'));
  });

// =======================
// Time settings
// =======================
var startDate = '2018-12-01';   // to include DJF 2019 (Dec 2018 – Feb 2019)
var endDate = '2025-05-31';     // until end of MAM 2025

// =======================
// NDVI & NDWI calculation
// =======================
function addIndices(image) {
  var ndvi = image.normalizedDifference(['B8', 'B4']).rename('NDVI');
  var ndwi = image.normalizedDifference(['B8', 'B11']).rename('NDWI');
  return image.addBands([ndvi, ndwi]);
}

// =======================
// Cloud masking (SCL-based)
// =======================
function maskS2clouds(image) {
  var scl = image.select('SCL');
  var mask = scl.neq(3).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));  // Exclude cloud and shadow
  return image.updateMask(mask);
}

// =======================
// S-2 image collection
// =======================
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterDate(startDate, endDate)
  .filterBounds(aoi)
  .map(maskS2clouds)
  .map(addIndices);

// =======================
// Annual median composites (2019–2024)
// =======================
var annualComposites = [];
for (var year = 2019; year <= 2024; year++) {
  var start = ee.Date.fromYMD(year, 1, 1);
  var end = ee.Date.fromYMD(year + 1, 1, 1);
  var yearIC = s2.filterDate(start, end);
  var composite = yearIC.median().clip(aoi).set({'year': year, 'season': 'annual'});
  if (composite.bandNames().size().getInfo() > 0) {
    annualComposites.push(composite);
  }
}
var annual_collection = ee.ImageCollection(annualComposites);

// =======================
// Seasonal median composites (2019 DJF to 2025 MAM)
// =======================
var seasons = [
  {name: 'DJF', months: [12, 1, 2]},
  {name: 'MAM', months: [3, 4, 5]},
  {name: 'JJA', months: [6, 7, 8]},
  {name: 'SON', months: [9, 10, 11]}
];

var seasonalComposites = [];
for (var year = 2019; year <= 2025; year++) {
  for (var i = 0; i < seasons.length; i++) {
    var season = seasons[i];
    var seasonName = season.name;
    var months = season.months;
    var start, end;

    if (seasonName === 'DJF') {
      start = ee.Date.fromYMD(year - 1, 12, 1);
      end   = ee.Date.fromYMD(year, 3, 1);
    } else if (year === 2025 && (seasonName === 'JJA' || seasonName === 'SON')) {
      continue; 
    } else {
      start = ee.Date.fromYMD(year, months[0], 1);
      var nextMonth = (months[2] % 12) + 1;
      end = ee.Date.fromYMD(year, nextMonth, 1);
    }

    var seasonIC = s2.filterDate(start, end);
    var composite = seasonIC.median().clip(aoi).set({'year': year, 'season': seasonName});
    if (composite.bandNames().size().getInfo() > 0) {
      seasonalComposites.push(composite);
    }
  }
}
var seasonal_collection = ee.ImageCollection(seasonalComposites);

// =======================
// Sample points from composites
// =======================
var fullBandList = [
  'B2', 'B3', 'B4', 'B5', 'B6', 'B7',
  'B8', 'B8A', 'B11', 'B12',
  'NDVI', 'NDWI'
];

var annualSamples = annual_collection.map(function(img) {
  var year = img.get('year');
  var season = img.get('season');
  return img.select(fullBandList)
    .sampleRegions({
      collection: all_points,
      properties: ['class', 'location_id'],
      scale: 20,
      tileScale: 4
    })
    .map(function(f) {
      return f.set({'year': year, 'season': season});
    });
}).flatten();

var seasonalSamples = seasonal_collection.map(function(img) {
  var year = img.get('year');
  var season = img.get('season');
  return img.select(fullBandList)
    .sampleRegions({
      collection: all_points,
      properties: ['class', 'location_id'],
      scale: 20,
      tileScale: 4
    })
    .map(function(f) {
      return f.set({'year': year, 'season': season});
    });
}).flatten();

print('Annual Sample size:', annualSamples.size());
print('Seasonal Sample size:', seasonalSamples.size());

// =======================
// 9. Export to Drive
// =======================
Export.table.toDrive({
  collection: annualSamples,
  description: 'test_thin43_s2_annual_bands2',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: seasonalSamples,
  description: 'test_thin43_s2_seasonal_bands2',
  fileFormat: 'CSV'
});
