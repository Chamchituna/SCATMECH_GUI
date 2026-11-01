import copy

from PyQt5.QtCore import Qt
from PyQt5.QtWidgets import QHBoxLayout, QMessageBox, QPushButton, QWidget


def node_expression(node):
    values = [f"{name}={value}" for name, value in node["params"].items()]
    return node["model"] + "(" + ",".join(values) + ")"


def prepare_profiles(raw, root_node):
    from custom_profile import compile_profile, profile_from_dict, profile_to_dict, remap_boundary_media
    from scatmech_gratings import CROSS_GRATING_SPECS, coerce_cross_grating, coerce_one_d_grating

    if not isinstance(raw, dict):
        raise ValueError("custom_profiles must be an object.")
    if not isinstance(root_node, dict):
        raise ValueError("The custom profile grating must be an object.")
    coerce = coerce_cross_grating if root_node.get("model") in CROSS_GRATING_SPECS else coerce_one_d_grating
    tree = coerce(root_node)
    documents = {}
    replacements = []
    targets = []
    for path, saved in sorted(raw.items(), key=lambda item: str(item[0])):
        if not isinstance(path, str) or not path.startswith("grating"):
            raise ValueError("Invalid custom profile target.")
        parts = path.split(".")
        if parts[0] != "grating" or (len(parts) - 1) % 2:
            raise ValueError(f"Invalid custom profile target: {path}")
        if any(path.startswith(previous + ".children.") for previous in documents):
            raise ValueError("Custom profile targets cannot overlap.")
        target = tree
        for index in range(1, len(parts), 2):
            if parts[index] != "children" or parts[index + 1] not in target.get("children", {}):
                raise ValueError(f"Custom profile target no longer exists: {path}")
            target = target["children"][parts[index + 1]]
        document = profile_from_dict(saved)
        compiled = compile_profile(document)
        if target.get("model") != compiled.node["model"]:
            raise ValueError(f"Custom profile does not match its grating model: {path}")
        replacements.append((target, compiled.node))
        targets.append((parts[2::2], compiled.node))
        documents[path] = profile_to_dict(document)
    tree = remap_boundary_media(tree, replacements)
    for path, node in targets:
        target = tree
        for name in path:
            target = target["children"][name]
        target.clear()
        target.update(copy.deepcopy(node))
    return documents, normalize_profile_boundaries(tree) if documents else tree


class CustomProfileBinding:
    def __init__(self, owner, dimension, get_node, set_node, context_provider, lock_fields, changed=None, media_provider=None):
        self.owner = owner
        self.dimension = dimension
        self.get_node = get_node
        self.set_node = set_node
        self.context_provider = context_provider
        self.lock_fields = lock_fields
        self.changed = changed
        self.media_provider = media_provider
        self.document = None
        self.editor = None
        self.updating = False
        self._compiled_node = None
        self.controls = QWidget(owner)
        layout = QHBoxLayout(self.controls)
        layout.setContentsMargins(0, 0, 0, 0)
        self.button = QPushButton("Custom Grating Profile…", self.controls)
        self.button.setToolTip(
            "Edit a 1D cross-section using material regions, surface points, or a native definition. New profiles start with a constrained core and spacers."
            if dimension == 1 else "Edit a 2D unit cell with circular, rectangular, or polygonal holes and posts."
        )
        self.detach_button = QPushButton("Use external file", self.controls)
        self.detach_button.setToolTip("Stop managing this geometry in the editor and unlock the native file parameters.")
        self.detach_button.hide()
        layout.addWidget(self.button)
        layout.addWidget(self.detach_button)
        self.button.clicked.connect(self._show_checked)
        self.detach_button.clicked.connect(self.detach)

    def _show_checked(self):
        try:
            self.show()
        except (ValueError, OSError) as exc:
            QMessageBox.warning(self.owner, "Custom Grating Profile", str(exc))

    def show(self):
        from custom_profile import default_editor_profile, profile_from_dict, profile_to_dict
        from custom_profile_editor import CustomProfileEditor

        if self.editor is not None and self.editor.isVisible():
            self.editor.raise_()
            self.editor.activateWindow()
            return self.editor
        params = self.media_provider() if self.media_provider is not None else self.get_node().get("params", {})
        if self.document is None:
            document = default_editor_profile(self.dimension, params.get("medium_i", "(1,0)"), params.get("medium_t", "(4.05,0.05)"))
        else:
            document = profile_from_dict(profile_to_dict(self.document))
            if self.media_provider is not None:
                from custom_profile import native_material_value
                from gui_validation import optical_constant

                for name in ("medium_i", "medium_t"):
                    inherited = params.get(name, getattr(document, name))
                    current = optical_constant(native_material_value(getattr(document, name)), name)
                    updated = optical_constant(native_material_value(inherited), name)
                    if updated != current:
                        setattr(document, name, inherited)
        if self.editor is None:
            self.editor = CustomProfileEditor(parent=self.owner.window(), document=document,
                context_provider=self.context_provider, apply_callback=self.apply)
            self.editor.setWindowModality(Qt.WindowModal)
            self.editor.setAttribute(Qt.WA_DeleteOnClose, False)
        else:
            self.editor.reload_document(document)
        self.editor.show()
        self.editor.raise_()
        self.editor.activateWindow()
        return self.editor

    def apply(self, document, node):
        from custom_profile import compile_profile, profile_from_dict, profile_to_dict

        if document.dimension != self.dimension:
            raise ValueError("The custom profile dimension does not match this grating.")
        saved = profile_to_dict(document)
        copied = profile_from_dict(saved)
        compiled = compile_profile(copied)
        if compiled.node != node:
            raise ValueError("The validated profile no longer matches the draft. Refresh the preview before applying.")
        previous = copy.deepcopy(self.get_node())
        self.updating = True
        try:
            try:
                self.set_node(copy.deepcopy(node))
            except Exception:
                self.set_node(previous)
                raise
            self.document = copied
            self._compiled_node = copy.deepcopy(node)
            self.update_controls()
        finally:
            self.updating = False
        if self.changed is not None:
            self.changed()

    def restore(self, saved, node):
        self.document = None
        if saved is not None:
            from custom_profile import profile_from_dict

            self.document = profile_from_dict(saved)
        self._compiled_node = copy.deepcopy(node) if self.document is not None else None
        self.update_controls()

    def detach(self, *_):
        if self.updating:
            return
        active = self.document is not None
        self.document = None
        self._compiled_node = None
        self.update_controls()
        if active and self.changed is not None:
            self.changed()

    def update_controls(self):
        active = self.document is not None
        self.lock_fields(active)
        self.detach_button.setVisible(active)

    def ensure(self):
        if self.document is None or self.updating:
            return
        from custom_profile import compile_profile

        self.updating = True
        try:
            compiled = compile_profile(self.document)
            self._compiled_node = copy.deepcopy(compiled.node)
            if self.get_node() != compiled.node:
                self.set_node(copy.deepcopy(compiled.node))
                self.update_controls()
        finally:
            self.updating = False

    def saved(self):
        if self.document is None:
            return None
        from custom_profile import profile_to_dict

        return profile_to_dict(self.document)

    def close(self):
        if self.editor is not None:
            self.editor.close()


def grating_editor_context(owner):
    context = {"wavelength_um": .532, "incident_deg": 0, "rotation_deg": 0, "configuration": 0}
    parent = owner.parentWidget()
    root = owner
    path = []
    while parent is not None:
        children = getattr(parent, "child_editors", {})
        name = next((name for name, child in children.items() if child is root), None)
        if name is not None:
            path.insert(0, name)
            root = parent
        provider = getattr(parent, "_custom_profile_context", None)
        if provider is not None:
            context.update(provider())
        parent = parent.parentWidget()
    if path:
        context["validation_grating"] = root.to_node()
        context["validation_path"] = path
        context["validation_profiles"] = editor_profiles(root)
    return context


def editor_profiles(editor, path="grating"):
    result = {}
    saved = editor.custom_profile_binding.saved()
    if saved is not None:
        result[path] = saved
    for name, child in getattr(editor, "child_editors", {}).items():
        result.update(editor_profiles(child, path + ".children." + name))
    return result


def restore_editor_profiles(editor, documents, path="grating"):
    editor.custom_profile_binding.restore(documents.get(path), editor.to_node())
    for name, child in getattr(editor, "child_editors", {}).items():
        restore_editor_profiles(child, documents, path + ".children." + name)


def grating_editor_media(owner):
    values = owner.to_node().get("params", {}).copy()
    parent = owner.parentWidget()
    while parent is not None:
        if getattr(parent, "cross_grating_editor", None) is owner:
            values["medium_i"] = "(1,0)"
            values["medium_t"] = parent.substrate.text()
            break
        parent = parent.parentWidget()
    return values


def has_editor_profiles(editor):
    return editor.custom_profile_binding.document is not None or any(
        has_editor_profiles(child) for child in getattr(editor, "child_editors", {}).values())


def normalize_profile_boundaries(node):
    from custom_profile import native_material_value

    result = copy.deepcopy(node)
    for name in ("medium_i", "medium_t"):
        if name in result.get("params", {}):
            result["params"][name] = native_material_value(result["params"][name])
    for name, child in result.get("children", {}).items():
        result["children"][name] = normalize_profile_boundaries(child)
    return result
