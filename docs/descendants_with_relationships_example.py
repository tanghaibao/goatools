#!/usr/bin/env python
"""
Example: Getting all descendants with relationships (part_of, regulates, etc.)

This example demonstrates how to get all descendants of a GO term including
those connected through optional relationships like 'part_of' and 'regulates'.

The functionality already exists in GOAtools - you just need to:
1. Load the GO DAG with optional_attrs={'relationship'}
2. Use the appropriate method depending on your needs
"""

from goatools.base import get_godag
from goatools.godag.go_tasks import get_go2descendants


def main():
    """Run the example demonstrating descendants with relationships."""
    # IMPORTANT: Load GO DAG with relationships
    # Use go.obo (not go-basic.obo) to include all relationships
    print("Loading GO DAG with relationships...")
    godag = get_godag('tests/data/heartjogging.obo', optional_attrs={'relationship'})

    # Example GO term
    GO_ID = 'GO:0003143'  # embryonic heart tube morphogenesis
    goterm = godag[GO_ID]

    print(f"\nExample GO term: {GO_ID} - {goterm.name}")
    print("=" * 70)

    # ============================================================================
    # METHOD 1: For a single GO term - using GOTerm methods
    # ============================================================================
    print("\n1. Using GOTerm methods (for single GO term):")
    print("-" * 70)

    # Get descendants through is_a only
    descendants_isa = goterm.get_all_children()
    print(f"   is_a only:              {sorted(descendants_isa)}")

    # Get descendants through is_a AND all relationships (part_of, regulates, etc.)
    descendants_all = goterm.get_all_lower()
    print(f"   is_a + all relations:   {sorted(descendants_all)}")

    # Show what we gained from relationships
    additional = descendants_all - descendants_isa
    print(f"   → Additional via relationships: {sorted(additional)}")

    # ============================================================================
    # METHOD 2: For multiple GO terms - using go_tasks functions
    # ============================================================================
    print("\n2. Using go_tasks functions (for batch processing):")
    print("-" * 70)

    # Get descendants for ALL GO terms in the DAG
    # Extract godag.values() once for efficiency
    goterms = godag.values()

    # Option A: is_a only
    go2descendants_isa = get_go2descendants(goterms, relationships=None)
    print(f"   is_a only:              {sorted(go2descendants_isa.get(GO_ID, set()))}")

    # Option B: is_a + ALL relationships
    go2descendants_all = get_go2descendants(goterms, relationships=True)
    print(f"   is_a + all relations:   {sorted(go2descendants_all.get(GO_ID, set()))}")

    # Option C: is_a + SELECTED relationships (e.g., only part_of)
    go2descendants_partof = get_go2descendants(goterms, relationships={'part_of'})
    print(f"   is_a + part_of only:    {sorted(go2descendants_partof.get(GO_ID, set()))}")

    # ============================================================================
    # METHOD 3: Using wr_hier.py command-line script
    # ============================================================================
    print("\n3. Using command-line script:")
    print("-" * 70)
    print("   # Without relationships (is_a only):")
    print("   $ goatools wr_hier GO:0003143 --dag=go.obo")
    print()
    print("   # With relationships (is_a + part_of, regulates, etc.):")
    print("   $ goatools wr_hier GO:0003143 --dag=go.obo -r")

    # ============================================================================
    # SUMMARY
    # ============================================================================
    print("\n" + "=" * 70)
    print("SUMMARY - How to get descendants with relationships:")
    print("=" * 70)
    print("""
1. MUST load GO DAG with: optional_attrs={'relationship'}
2. MUST use go.obo (not go-basic.obo) to include all relationships

For single GO term:
  - goterm.get_all_children()  -> descendants via is_a only
  - goterm.get_all_lower()     -> descendants via is_a + all relationships

For multiple GO terms:
  - get_go2descendants(terms, relationships=None)         -> is_a only
  - get_go2descendants(terms, relationships=True)         -> is_a + all
  - get_go2descendants(terms, relationships={'part_of'})  -> is_a + selected

For visualization:
  - goatools wr_hier GO_ID --dag=go.obo -r  -> hierarchy with relationships

See notebooks/children_and_descendants.ipynb for more examples!
""")


if __name__ == "__main__":
    main()
