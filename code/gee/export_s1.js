// ==========================
// AOI & Points
// ==========================
var aoi = ee.FeatureCollection('projects/ee-yudong9607/assets/paper/Buvuma_down2018');
var all_points = ee.FeatureCollection('projects/ee-yudong9607/assets/paper/Test_thin43_2')
  .map(function(f) { return f.set('location_id', f.get('down_id')); });

// ==========================
// Refined Lee + Ratio in dB
// ==========================
function toNatural(imgDB) {
  return ee.Image(10).pow(ee.Image(imgDB).divide(10));
}
function toDB(imgLin) {
  return ee.Image(imgLin).log10().multiply(10);
}
function RefinedLee(imgDB) {
  var myimg = toNatural(imgDB);
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3,3, weights3, 1, 1, false);
  var mean3 = myimg.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance3 = myimg.reduceNeighborhood(ee.Reducer.variance(), kernel3);

  var sample_weights = ee.List([
    [0,0,0,0,0,0,0],[0,1,0,1,0,1,0],[0,0,0,0,0,0,0],
    [0,1,0,1,0,1,0],[0,0,0,0,0,0,0],[0,1,0,1,0,1,0],[0,0,0,0,0,0,0]
  ]);
  var sample_kernel = ee.Kernel.fixed(7,7, sample_weights, 3,3, false);
  var sample_mean = mean3.neighborhoodToBands(sample_kernel);
  var sample_var  = variance3.neighborhoodToBands(sample_kernel);
  var gradients = sample_mean.select(1).subtract(sample_mean.select(7)).abs()
    .addBands(sample_mean.select(6).subtract(sample_mean.select(2)).abs())
    .addBands(sample_mean.select(3).subtract(sample_mean.select(5)).abs())
    .addBands(sample_mean.select(0).subtract(sample_mean.select(8)).abs());
  var max_gradient = gradients.reduce(ee.Reducer.max());
  var gradmask = gradients.eq(max_gradient).addBands(gradients.eq(max_gradient));
  var directions = sample_mean.select(1).subtract(sample_mean.select(4))
    .gt(sample_mean.select(4).subtract(sample_mean.select(7))).multiply(1)
    .addBands(sample_mean.select(6).subtract(sample_mean.select(4))
    .gt(sample_mean.select(4).subtract(sample_mean.select(2))).multiply(2))
    .addBands(sample_mean.select(3).subtract(sample_mean.select(4))
    .gt(sample_mean.select(4).subtract(sample_mean.select(5))).multiply(3))
    .addBands(sample_mean.select(0).subtract(sample_mean.select(4))
    .gt(sample_mean.select(4).subtract(sample_mean.select(8))).multiply(4))
    .addBands(sample_mean.select(0).multiply(0).add(5))
    .addBands(sample_mean.select(0).multiply(0).add(6))
    .addBands(sample_mean.select(0).multiply(0).add(7))
    .addBands(sample_mean.select(0).multiply(0).add(8));
  directions = directions.updateMask(gradmask).reduce(ee.Reducer.sum());

  var sample_stats = sample_var.divide(sample_mean.multiply(sample_mean));
  var sigmaV = sample_stats.toArray().arraySort().arraySlice(0,0,5).arrayReduce(ee.Reducer.mean(), [0]);

  var rect_weights = ee.List.repeat(ee.List.repeat(0,7),3).cat(ee.List.repeat(ee.List.repeat(1,7),4));
  var diag_weights = ee.List([
    [1,0,0,0,0,0,0],[1,1,0,0,0,0,0],[1,1,1,0,0,0,0],
    [1,1,1,1,0,0,0],[1,1,1,1,1,0,0],[1,1,1,1,1,1,0],[1,1,1,1,1,1,1]
  ]);
  var rect_kernel = ee.Kernel.fixed(7,7, rect_weights, 3, 3, false);
  var diag_kernel = ee.Kernel.fixed(7,7, diag_weights, 3, 3, false);
  var dir_mean = myimg.reduceNeighborhood(ee.Reducer.mean(), rect_kernel).updateMask(directions.eq(1));
  var dir_var  = myimg.reduceNeighborhood(ee.Reducer.variance(), rect_kernel).updateMask(directions.eq(1));
  dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), diag_kernel).updateMask(directions.eq(2)));
  dir_var  = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), diag_kernel).updateMask(directions.eq(2)));
  for (var i=1; i<4; i++) {
    dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_var  = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), rect_kernel.rotate(i)).updateMask(directions.eq(2*i+1)));
    dir_mean = dir_mean.addBands(myimg.reduceNeighborhood(ee.Reducer.mean(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
    dir_var  = dir_var.addBands(myimg.reduceNeighborhood(ee.Reducer.variance(), diag_kernel.rotate(i)).updateMask(directions.eq(2*i+2)));
  }
  dir_mean = dir_mean.reduce(ee.Reducer.sum());
  dir_var  = dir_var.reduce(ee.Reducer.sum());
  var varX = dir_var.subtract(dir_mean.multiply(dir_mean).multiply(sigmaV)).divide(sigmaV.add(1.0));
  var b = varX.divide(dir_var);
  var result = dir_mean.add(b.multiply(myimg.subtract(dir_mean)));
  return toDB(result.arrayGet(0)).rename('filter');
}
function refinedLeeBand(img, band, outName) {
  var filtered = RefinedLee(img.select([band]));
  return img.addBands(filtered.rename(outName), null, true);
}
function addFilteredBandsAndRatio(img) {
  var out = refinedLeeBand(img, 'VV', 'VV');
  out = refinedLeeBand(out, 'VH', 'VH');
  var vvvh = out.select('VV').subtract(out.select('VH')).rename('VVVH');
  return out.addBands(vvvh);
}

// ==========================
// Time settings & seasons
// ==========================
var startYear = 2019;
var endYear   = 2025;

var seasons = [
  {name: 'DJF', months: [12,1,2]},
  {name: 'MAM', months: [3,4,5]},
  {name: 'JJA', months: [6,7,8]},
  {name: 'SON', months: [9,10,11]}
];

// ==========================
// S1 collection
// ==========================
function s1Collection(start, end) {
  return ee.ImageCollection('COPERNICUS/S1_GRD')
    .filterBounds(aoi)
    .filterDate(start, end)
    .filter(ee.Filter.eq('instrumentMode', 'IW'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VV'))
    .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
    .select(['VV','VH'])
    .map(function(img) {
        return img.set('orbit_pass', img.get('orbitProperties_pass'));
    })
    .map(addFilteredBandsAndRatio);
}

// ==========================
// Annual composites (2019–2024)
// ==========================
var annualComposites = [];
for (var y = 2019; y <= 2024; y++) {
  var start = ee.Date.fromYMD(y, 1, 1);
  var end   = ee.Date.fromYMD(y + 1, 1, 1);
  var s1 = s1Collection(start, end);
  var comp = s1.median().clip(aoi).set({'year': y, 'season': 'annual'});
  if (comp.bandNames().size().getInfo() > 0) {
    annualComposites.push(comp);
  }
}
var annual_collection = ee.ImageCollection(annualComposites);

// ==========================
// Seasonal composites (DJF 2019 → MAM 2025)
// ==========================
var seasonalComposites = [];
for (var year = startYear; year <= endYear; year++) {
  for (var i = 0; i < seasons.length; i++) {
    var s = seasons[i];
    var seasonName = s.name;
    var months = s.months;
    var start, end;
    if (seasonName === 'DJF') {
      start = ee.Date.fromYMD(year - 1, 12, 1);
      end   = ee.Date.fromYMD(year, 3, 1);
    } else if (year === 2025 && (seasonName === 'JJA' || seasonName === 'SON')) {
      continue;  // Only process up to MAM 2025
    } else {
      start = ee.Date.fromYMD(year, months[0], 1);
      var nextMonth = (months[2] % 12) + 1;
      end = ee.Date.fromYMD(year, nextMonth, 1);
    }

    var s1s = s1Collection(start, end);
    var comp = s1s.median().clip(aoi).set({'year': year, 'season': seasonName});
    if (comp.bandNames().size().getInfo() > 0) {
      seasonalComposites.push(comp);
    }
  }
}
var seasonal_collection = ee.ImageCollection(seasonalComposites);

// ==========================
// Sample & Export
// ==========================
var bandList = ['VV','VH','VVVH'];

var annualSamples = annual_collection.map(function(img) {
  var year = img.get('year');
  return img.select(bandList)
    .sampleRegions({
      collection: all_points,
      properties: ['class', 'location_id'],
      scale: 10,
      tileScale: 4
    }).map(function(f) {
      return f.set({'year': year, 'season': 'annual'});
    });
}).flatten();

var seasonalSamples = seasonal_collection.map(function(img) {
  var year = img.get('year');
  var season = img.get('season');
  return img.select(bandList)
    .sampleRegions({
      collection: all_points,
      properties: ['class', 'location_id'],
      scale: 10,
      tileScale: 4
    }).map(function(f) {
      return f.set({'year': year, 'season': season});
    });
}).flatten();

Export.table.toDrive({
  collection: annualSamples,
  description: 'test_thin43_s1_annual_filtered2',
  fileFormat: 'CSV'
});

Export.table.toDrive({
  collection: seasonalSamples,
  description: 'test_thin43_s1_seasonal_filtered2',
  fileFormat: 'CSV'
});
