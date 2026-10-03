"""The Tk dialog: create it and re-lay it out for every choice that drives the
layout, checking what is shown. Skipped when no display is available."""

import pytest


@pytest.fixture()
def tk_node():
    import tkinter as tk

    try:
        root = tk.Tk()
    except tk.TclError:
        pytest.skip("no display available for Tk")
    root.withdraw()
    import Pmw
    import seamm

    Pmw.initialise(root)
    flowchart = seamm.Flowchart(namespace="org.molssi.seamm", directory=".")
    tk_flowchart = seamm.TkFlowchart(
        master=root, flowchart=flowchart, namespace="org.molssi.seamm.tk"
    )
    node = flowchart.create_node("MBE")
    flowchart.add_node(node)
    plugin = tk_flowchart.plugin_manager.get("MBE")
    tk_node = plugin.create_tk_node(
        tk_flowchart=tk_flowchart, node=node, canvas=tk_flowchart.canvas, x=100, y=100
    )
    yield tk_node
    root.destroy()


def mapped(tk_node, key):
    return tk_node[key].grid_info() != {}


def test_dialog_layouts(tk_node):
    tk_node.create_dialog()
    tk_node.reset_dialog()
    # Defaults: order 3, single cutoffs, connected triples, no periodic level
    assert mapped(tk_node, "pair cutoff") and mapped(tk_node, "triple cutoff")
    assert not mapped(tk_node, "pair cutoff table")
    assert not mapped(tk_node, "periodic monomers")
    assert not mapped(tk_node, "periodic ranks")

    # The order caps at 3 in the GUI
    assert list(tk_node["maximum order"].combobox.cget("values")) == ["1", "2", "3"]

    tk_node["maximum order"].set("2")
    tk_node.reset_dialog()
    assert mapped(tk_node, "pair cutoff")
    assert not mapped(tk_node, "triple rule")
    assert not mapped(tk_node, "triple cutoff")

    tk_node["maximum order"].set("1")
    tk_node.reset_dialog()
    assert not mapped(tk_node, "cutoffs")
    assert not mapped(tk_node, "pair cutoff")

    tk_node["maximum order"].set("3")
    tk_node["triple rule"].set("none")
    tk_node.reset_dialog()
    assert mapped(tk_node, "triple rule") and not mapped(tk_node, "triple cutoff")

    tk_node["triple rule"].set("compact")
    tk_node["cutoffs"].set("table by type pair")
    tk_node.reset_dialog()
    assert mapped(tk_node, "pair cutoff table")
    assert mapped(tk_node, "triple cutoff table")
    assert not mapped(tk_node, "pair cutoff")
    assert not mapped(tk_node, "triple cutoff")

    tk_node["periodic low level"].set("XNN:MLFF@something")
    tk_node.reset_dialog()
    assert mapped(tk_node, "periodic monomers")
    assert mapped(tk_node, "periodic pair cutoff")
    assert mapped(tk_node, "periodic ranks")

    tk_node["periodic low level"].set("none")
    tk_node.reset_dialog()
    assert not mapped(tk_node, "periodic pair cutoff")

    for systems in ("all", "current"):
        tk_node["source systems"].set(systems)
        tk_node.reset_dialog()
    tk_node["source configurations"].set("name is")
    tk_node.reset_dialog()
    assert mapped(tk_node, "source configuration name")


def test_counterpoise_is_not_offered(tk_node):
    """Counterpoise arrives in phase 3; it must not appear at all yet."""
    assert not any("counterpoise" in key for key in tk_node.node.parameters)


def test_variable_order_does_not_break_the_layout(tk_node):
    """A typed $variable in the order field must not make reset_dialog raise."""
    tk_node.create_dialog()
    tk_node["maximum order"].set("$order")
    tk_node.reset_dialog()
    assert tk_node["triple rule"].grid_info() != {}
