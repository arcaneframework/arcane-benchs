#!/bin/bash

export WORKSPACE_DIR="/workspaces"
export FOLDER_DIR="$WORKSPACE_DIR/arcane-benchs"

if [[ -e "$FOLDER_DIR/.vscode/c_cpp_properties.json" ]] && [[ ! -e "$FOLDER_DIR/.vscode/c_cpp_properties_host.json" ]]
then
  mv "$FOLDER_DIR/.vscode/c_cpp_properties.json" "$FOLDER_DIR/.vscode/c_cpp_properties_host.json"
fi

cp "$FOLDER_DIR/.devcontainer/c_cpp_properties.json" "$FOLDER_DIR/.vscode/c_cpp_properties.json"
