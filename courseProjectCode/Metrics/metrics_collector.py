

import ast
import json
import os
import re
from typing import Dict, List, Tuple


ROOT_DIR = os.getcwd()

SKIP_DIRS = {
    "node_modules",
    "courseProjectDocs",
    "courseProjectCode", 
    ".git",
    "__pycache__",
}

SOURCE_EXTENSIONS = {".py"}


def count_python_functions(file_content: str) -> Tuple[int, List[Tuple[int, int]]]:
    try:
        tree = ast.parse(file_content)
    except SyntaxError:
        return 0, []

    function_spans = []
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            # end_lineno is available in Python 3.8+
            start_line = getattr(node, "lineno", None)
            end_line = getattr(node, "end_lineno", None)
            if start_line is not None and end_line is not None:
                function_spans.append((start_line, end_line))
    return len(function_spans), function_spans


def count_js_functions(file_content: str) -> Tuple[int, List[Tuple[int, int]]]:
    lines = file_content.splitlines()
    count = 0
    spans = []
    for idx, line in enumerate(lines, start=1):
        stripped = line.strip()
        if stripped.startswith("//") or stripped.startswith("/*"):
            continue
        if re.search(r"\bfunction\b", stripped) or re.search(r"=>", stripped):
            count += 1
            spans.append((idx, idx))
    return count, spans


def approximate_cyclomatic_complexity(lines: List[str]) -> int:
    complexity = 1  # Base complexity
    decision_keywords = [
        "if ", "for ", "while ", "case ", "switch ", "catch ", "&&", "||", "?",
        "elif ", "except ",
    ]
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#") or stripped.startswith("//"):
            continue
        for keyword in decision_keywords:
            if keyword in stripped:
                complexity += 1
                break
    return complexity


def analyse_file(filepath: str) -> Dict[str, object]:
    try:
        with open(filepath, "r", encoding="utf-8", errors="ignore") as f:
            content = f.read()
    except (OSError, UnicodeDecodeError):
        return {}

    lines = content.splitlines()
    code_lines = 0
    comment_lines = 0
    in_block_comment = False

    for line in lines:
        stripped = line.strip()
        if not stripped:
            continue
        if in_block_comment:
            comment_lines += 1
            if "*/" in stripped:
                in_block_comment = False
            continue
        if stripped.startswith("/*"):
            comment_lines += 1
            if "*/" not in stripped:
                in_block_comment = True
            continue
        if stripped.startswith("#") or stripped.startswith("//"):
            comment_lines += 1
            continue
        if stripped.startswith("\"\"\""):
            comment_lines += 1
            continue
        code_lines += 1

    ext = os.path.splitext(filepath)[1]
    functions_count = 0
    function_spans: List[Tuple[int, int]] = []
    if ext == ".py":
        functions_count, function_spans = count_python_functions(content)
    elif ext == ".js":
        functions_count, function_spans = count_js_functions(content)

    total_function_lines = 0
    for start, end in function_spans:
        if end >= start:
            total_function_lines += end - start + 1
    average_function_length = (
        (total_function_lines / functions_count) if functions_count > 0 else 0
    )

    complexity = approximate_cyclomatic_complexity(lines)

    parts = filepath.lower().split(os.sep)
    is_test_file = any(
        part.startswith("test") for part in parts if part not in {"", "."}
    )

    test_functions_count = 0
    if is_test_file:
        if ext == ".py":
            try:
                tree = ast.parse(content)
            except SyntaxError:
                tree = None
            if tree is not None:
                for node in ast.walk(tree):
                    if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                        if node.name.startswith("test"):
                            test_functions_count += 1
        elif ext == ".js":
            test_functions_count = len(re.findall(r"\b(it|describe)\s*\(", content))

    return {
        "file": filepath,
        "lines_of_code": code_lines,
        "comment_lines": comment_lines,
        "comment_ratio": (comment_lines / code_lines) if code_lines > 0 else 0,
        "functions": functions_count,
        "average_function_length": average_function_length,
        "cyclomatic_complexity": complexity,
        "is_test_file": is_test_file,
        "test_functions": test_functions_count,
    }


def walk_repository(root_dir: str) -> List[Dict[str, object]]:
    results = []
    for dirpath, dirnames, filenames in os.walk(root_dir):
        # Remove skipped directories from traversal
        dirnames[:] = [d for d in dirnames if d not in SKIP_DIRS]
        for filename in filenames:
            ext = os.path.splitext(filename)[1]
            if ext in SOURCE_EXTENSIONS:
                filepath = os.path.join(dirpath, filename)
                metrics = analyse_file(filepath)
                if metrics:
                    results.append(metrics)
    return results


def aggregate_metrics(results: List[Dict[str, object]]) -> Dict[str, object]:

    total_code_lines = sum(item["lines_of_code"] for item in results)
    total_comment_lines = sum(item["comment_lines"] for item in results)
    total_functions = sum(item["functions"] for item in results)
    total_complexity = sum(item["cyclomatic_complexity"] for item in results)
    total_files = len(results)
    
    total_function_lines = sum(
        item["average_function_length"] * item["functions"] for item in results
    )
    average_function_length = (
        total_function_lines / total_functions if total_functions > 0 else 0
    )
    comment_ratio = (
        (total_comment_lines / total_code_lines) if total_code_lines > 0 else 0
    )

    test_files = [item for item in results if item["is_test_file"]]
    total_test_files = len(test_files)
    total_test_lines = sum(item["lines_of_code"] for item in test_files)
    total_test_functions = sum(item["test_functions"] for item in test_files)
    test_ratio = (
        (total_test_lines / total_code_lines) if total_code_lines > 0 else 0
    )

    aggregated = {
        "total_files": total_files,
        "total_code_lines": total_code_lines,
        "total_comment_lines": total_comment_lines,
        "comment_ratio": comment_ratio,
        "total_functions": total_functions,
        "average_function_length": average_function_length,
        "total_cyclomatic_complexity": total_complexity,
        "total_test_files": total_test_files,
        "total_test_lines": total_test_lines,
        "total_test_functions": total_test_functions,
        "test_ratio": test_ratio,
    }
    return aggregated


def main() -> None:
    results = walk_repository(ROOT_DIR)
    aggregated = aggregate_metrics(results)
    report = {
        "files": results,
        "summary": aggregated,
    }
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()