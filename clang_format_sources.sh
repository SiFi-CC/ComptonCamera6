find . \( -iname "*.hh" -or -iname "*.cc" -or -iname "*.C" \) -exec clang-format --style=file -i {} \;
