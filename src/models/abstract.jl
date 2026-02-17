"""
Abstract supertype for barcode models.

Subtypes must expose a `params::ModelParams` field; `simulate` relies on this.
"""
abstract type AbstractBarcodeModel end
