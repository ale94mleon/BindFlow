import subprocess
import difflib
import os

# -------------------------
# CONFIG: Set your commit
# -------------------------
OLD_COMMIT = "33fc6d3a0e93ffa74852f2721550de6249ac7d2e"  # replace with your old commit
# -------------------------

def get_files_at_commit(commit):
    """Get a list of all files in a given commit."""
    result = subprocess.run(
        ["git", "ls-tree", "-r", "--name-only", commit],
        capture_output=True, text=True
    )
    return result.stdout.strip().split("\n")

def get_file_content(commit, filepath):
    """Get the content of a file at a specific commit."""
    result = subprocess.run(
        ["git", "show", f"{commit}:{filepath}"],
        capture_output=True
    )
    if result.returncode != 0:
        return []

    try:
        text = result.stdout.decode("utf-8")
        return text.splitlines()
    except UnicodeDecodeError:
        return []  # skip binary files

def get_current_file_content(filepath):
    """Get the content of a file in the working tree."""
    if not os.path.exists(filepath):
        return []
    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
        return f.read().splitlines()

def similarity_ratio(lines_old, lines_new):
    """Compute similarity ratio using difflib."""
    sm = difflib.SequenceMatcher(None, lines_old, lines_new)
    return sm.ratio()  # 0.0 -> totally different, 1.0 -> identical

def main():
    old_files = get_files_at_commit(OLD_COMMIT)
    
    # get list of all current files
    current_files = []
    for root, dirs, files in os.walk("."):
        for f in files:
            current_files.append(os.path.relpath(os.path.join(root, f), "."))

    total_old_lines = 0
    weighted_similarity = 0

    for old_file in old_files:
        old_lines = get_file_content(OLD_COMMIT, old_file)
        if not old_lines:
            continue  # skip binary or missing file

        total_old_lines += len(old_lines)

        # find best matching current file
        best_ratio = 0
        for current_file in current_files:
            new_lines = get_current_file_content(current_file)
            if not new_lines:
                continue
            ratio = similarity_ratio(old_lines, new_lines)
            if ratio > best_ratio:
                best_ratio = ratio

        # weight by number of lines in the old file
        weighted_similarity += best_ratio * len(old_lines)

    if total_old_lines == 0:
        print("No text files found in the old commit!")
        return

    percentage_remaining = (weighted_similarity / total_old_lines) * 100
    print(f"Approximate percentage of original code remaining: {percentage_remaining:.2f}%")

if __name__ == "__main__":
    main()