var aoi = ee.FeatureCollection('projects/ee-yudong9607/assets/paper/Buvuma_down2018');

var all_points = ee.FeatureCollection('projects/ee-yudong9607/assets/paper/Test_thin43_2')
  .map(function(f) {
    return f.set('location_id', f.get('down_id'));
  });

// ==========================
// Load Copernicus GLO-30 DEM (mosaic) and compute slope/aspect
// ==========================
var dem_raw = ee.ImageCollection("COPERNICUS/DEM/GLO30").mosaic();
var dem = dem_raw.select('DEM');  // Only the elevation band
var slope = ee.Terrain.slope(dem);
var aspect = ee.Terrain.aspect(dem);

var topo = dem.rename('elevation')
  .addBands(slope.rename('slope'))
  .addBands(aspect.rename('aspect'));


var topo = dem.rename('elevation')
  .addBands(slope.rename('slope'))
  .addBands(aspect.rename('aspect'));

// ==========================
// Sample topographic values at points
// ==========================
var topoSamples = topo.sampleRegions({
  collection: all_points,
  properties: ['location_id', 'class'],
  scale: 30,
  tileScale: 4
});

// ==========================
// Export
// ==========================
Export.table.toDrive({
  collection: topoSamples,
  description: 'test_thin43_topo2',
  fileFormat: 'CSV'
});
