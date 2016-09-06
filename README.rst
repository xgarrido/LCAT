Falaise (Light)CAT plugin
=========================

This plugin for Falaise implement some processing modules for the SuperNEMO
reconstruction pipeline.

It provides a light version of the CellularAutomatonTracker (CAT) written
by F. Nova, a library dedicated to the clustering of tracker Geiger hits in the
SuperNEMO experiment. This version only focuses on clustering and all
histogramming, track fitting and other physical interpretation are removed.

Unlike the original work, the LCAT library requires Bayeux library.
