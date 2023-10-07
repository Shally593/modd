import { Component } from '@angular/core';
import * as L from 'leaflet';

@Component({
  selector: 'app-home',
  templateUrl: 'home.page.html',
  styleUrls: ['home.page.scss'],
})
export class HomePage {
  map!: L.Map;

  constructor() {}
  // ngOnInit() {

  // }

  ionViewDidEnter() {
    this.map = L.map('mapId').setView([-7.781497658743231, 110.36722982606936], 10)

    L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
      attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
    }).addTo(this.map);

    const baseMaps = {
      'OpenStreetMap' : L.tileLayer('https://tile.openstreetmap.org/{z}/{x}/{y}.png', {
        attribution: '&copy; <a href="https://www.openstreetmap.org/copyright">OpenStreetMap</a> contributors'
      }),

      'Stadia' : L.tileLayer('https://server.arcgisonline.com/ArcGIS/rest/services/World_Imagery/MapServer/tile/{z}/{y}/{x}', {
        attribution: 'Tiles &copy; Esri &mdash; Source: Esri, i-cubed, USDA, USGS, AEX, GeoEye, Getmapping, Aerogrid, IGN, IGP, UPR-EGP, and the GIS User Community'
      })
         
    };
    L.control.layers(baseMaps).addTo(this.map);


    // L.marker([-7.781497658743231, 110.36722982606936]).addTo(this.map)
    // .bindPopup('Marker1')
    // .openPopup();

    // L.marker([-7.723262973912956, 110.38033435824343]).addTo(this.map)
    // .bindPopup('Marker1')
    // .openPopup();

    const markerIcon = L.icon({
      iconUrl: 'https://unpkg.com/leaflet@1.7.1/dist/images/marker-icon.png', // Ganti dengan URL ikon marker default dari CDN
      iconRetinaUrl: 'https://unpkg.com/leaflet@1.7.1/dist/images/marker-icon-2x.png', // Ganti dengan URL ikon marker default 2x dari CDN
      shadowUrl: 'https://unpkg.com/leaflet@1.7.1/dist/images/marker-shadow.png', // Ganti dengan URL bayangan marker default dari CDN
      iconSize: [25, 41], // Sesuaikan dengan ukuran ikon Anda
      iconAnchor: [12, 41], // Sesuaikan dengan titik penunjuk ikon Anda
    });
    
    
    const marker = L.marker([-7.723262973912956, 110.38033435824343], { icon: markerIcon }).addTo(this.map).bindPopup('Marker1').openPopup();
    const marker1 = L.marker([-7.781497658743231, 110.36722982606936], { icon: markerIcon }).addTo(this.map).bindPopup('Marker2').openPopup();
  }

    // Add layers to map
    // streetMapLayer.addTo(this.map);

    // // Layer control
    // const baseMaps = {
    //   'Street Map': streetMapLayer,
    //   'Satellite Map': satelliteMapLayer,
    //   'Terrain Map': terrainMapLayer,
    //   'Watercolor Map': watercolorMapLayer,
    //   'Dark Map': darkMapLayer
    // };

    // L.control.layers(baseMaps).addTo(this.map);
  }