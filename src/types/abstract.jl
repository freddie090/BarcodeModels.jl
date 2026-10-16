"""
Abstract supertype for barcode models.

Subtypes must expose a `params` field whose concrete type is model-specific.
"""
abstract type AbstractBarcodeModel end

"""Hybrid ODE/jump model family."""
abstract type HybridModel <: AbstractBarcodeModel end

"""Agent-based model family."""
abstract type ABMModel <: AbstractBarcodeModel end


