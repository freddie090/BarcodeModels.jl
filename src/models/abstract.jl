"""
Abstract supertype for barcode models.

Subtypes must expose a `params::ModelParams` field.
"""
abstract type AbstractBarcodeModel end

"""Hybrid ODE/jump model family."""
abstract type HybridModel <: AbstractBarcodeModel end

"""Agent-based model family."""
abstract type ABMModel <: AbstractBarcodeModel end
