//Copyright (c) 2016 eXact Lab s.r.l.
//Author: Gianfranco Gallizia <gianfranco.gallizia@exact-lab.it>

//Poly selector initialization

//Use Stamen tile server for general features
var layer = new L.StamenTileLayer("toner-lite");

//Create Bathymetry layer
var blayer = L.tileLayer(' http://maps.ngdc.noaa.gov/arcgis/rest/services/web_mercator/gebco08_contours/MapServer/tile/{z}/{y}/{x}', { maxNativeZoom: 8 });

//Create Leaflet map object (centered on the Mediterranean Sea)
var map = L.map('map').setView([39.00, 18.00], 5);

//Add Tile layers to map
map.addLayer(layer);

//Create Layer control
var layersKnob = L.control.layers({"Base": layer}, {"Bathymetry": blayer});

//Add Layer control to the map
layersKnob.addTo(map);

//Create FeatureGroup to store the drawn polygons
var drawnItems = new L.FeatureGroup();
map.addLayer(drawnItems);

// Initialize the draw control and pass it the FeatureGroup of editable layers
var drawControl = new L.Control.Draw({
    draw: { polyline: false, rectangle: false, circle: false, marker: false },
    edit: { featureGroup: drawnItems }
});

function refreshPointsList(points) {
    $('#points').empty();
    for (i in points) {
        $("#points").append("" + points[i].lat + " " + points[i].lng + "<br/>\n");
    }
    if (typeof exportPointsList != 'undefined')
        exportPointsList(points);
}

map.addControl(drawControl);
map.on('draw:created', function(e) {
    var layer = e.layer;

    drawnItems.addLayer(layer);
    refreshPointsList(layer.getLatLngs());
});

map.on('draw:edited', function(e) {
    var layers = e.layers.getLayers();
    //Take only the last one
    var layer = layers[layers.length-1];

    refreshPointsList(layer.getLatLngs());
});

map.on('draw:deleted', function(e) {
    $('#points').empty();
    $('#export').empty();
    $('#ui').empty();
});
