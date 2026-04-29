#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 1 ]]; then
  cat <<'EOF'
Usage: src/scripts/run_utilities.sh [--rebuild-image] [--no-build-cache] <command> [command args...]

Commands:
  compare-columns <input-file> [compare options]
  filter-columns  <input-file> <output-file> [filter options]
  filter-rows     <input-file> <output-file> [row filter options]
  replace-rows    <output-file> <input-file>... [replace options]
  merge-columns   <base-file> <extra-file> <output-file> [merge options]
  reorder-columns <input-file> <output-file> --column-order <cols>

Examples:
  src/scripts/run_utilities.sh compare-columns in.tsv --col-a mapped_hgvs_g --col-b ref_hgvs_g --output diffs.tsv
  src/scripts/run_utilities.sh filter-columns in.tsv out.tsv --keep-col a --keep-col b
  src/scripts/run_utilities.sh filter-rows in.tsv out.tsv --value-col a,b --match any
  src/scripts/run_utilities.sh replace-rows replaced.tsv base.tsv patch.tsv --key-col id
  src/scripts/run_utilities.sh merge-columns base.tsv extra.tsv merged.tsv --key-col id --add-col score
  src/scripts/run_utilities.sh reorder-columns in.tsv out.tsv --column-order id,gene,value

Notes:
  - Paths are interpreted relative to /work in the container.
  - By default /work maps to ./data on the host.
  - Override mount directory with VARIANT_DATA_DIR=/absolute/path.
  - Images are reused by default for fast runs.
  - Add --rebuild-image to force rebuilding the tools image.
  - Add --no-build-cache with --rebuild-image for a clean rebuild.
EOF
  exit 1
fi

compose_build_flag=""
compose_no_cache_flag=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --rebuild-image)
      compose_build_flag="--build"
      shift
      ;;
    --no-build-cache)
      compose_no_cache_flag="--no-cache"
      shift
      ;;
    --)
      shift
      break
      ;;
    *)
      break
      ;;
  esac
done

if [[ $# -lt 1 ]]; then
  echo "error: missing command (expected: filter-columns, filter-rows, replace-rows, merge-columns, reorder-columns, or compare-columns)" >&2
  exit 2
fi

command_name="$1"
shift

map_to_container_path() {
  local path="$1"
  local prefer_app_path="${2:-auto}"

  if [[ "$path" == /work/* || "$path" == /usr/src/app/* ]]; then
    printf '%s\n' "$path"
    return
  fi

  if [[ "$prefer_app_path" == "app" ]]; then
    printf '/usr/src/app/%s\n' "$path"
    return
  fi

  if [[ "$prefer_app_path" == "work" ]]; then
    printf '/work/%s\n' "$path"
    return
  fi

  if [[ -f "$path" ]]; then
    printf '/usr/src/app/%s\n' "$path"
  else
    printf '/work/%s\n' "$path"
  fi
}

mapped_args=()

case "$command_name" in
  compare-columns)
    if [[ $# -lt 1 ]]; then
      echo "error: $command_name requires <input-file> [compare options]" >&2
      exit 2
    fi
    input_path="$1"
    shift

    input_in_container="$(map_to_container_path "$input_path")"

    mapped_args=("$command_name" "$input_in_container")
    if [[ $# -gt 0 ]]; then
      mapped_args+=("$@")
    fi
    ;;

  filter-columns|filter-rows)
    if [[ $# -lt 2 ]]; then
      echo "error: $command_name requires <input-file> <output-file>" >&2
      exit 2
    fi
    input_path="$1"
    output_path="$2"
    shift 2

    input_in_container="$(map_to_container_path "$input_path")"
    if [[ "$output_path" == /work/* || "$output_path" == /usr/src/app/* ]]; then
      output_in_container="$output_path"
    elif [[ "$input_in_container" == /usr/src/app/* ]]; then
      output_in_container="/usr/src/app/$output_path"
    else
      output_in_container="/work/$output_path"
    fi

    mapped_args+=("$command_name" "$input_in_container" "$output_in_container")
    if [[ $# -gt 0 ]]; then
      mapped_args+=("$@")
    fi
    ;;

  replace-rows)
    if [[ $# -lt 2 ]]; then
      echo "error: replace-rows requires <output-file> <input-file>..." >&2
      exit 2
    fi

    output_path="$1"
    shift

    input_paths=()
    while [[ $# -gt 0 ]]; do
      if [[ "$1" == --* ]]; then
        break
      fi
      input_paths+=("$1")
      shift
    done

    if [[ ${#input_paths[@]} -lt 1 ]]; then
      echo "error: replace-rows requires at least one input file" >&2
      exit 2
    fi

    mapped_inputs=()
    for in_path in "${input_paths[@]}"; do
      mapped_inputs+=("$(map_to_container_path "$in_path")")
    done

    first_input="${mapped_inputs[0]}"
    if [[ "$output_path" == /work/* || "$output_path" == /usr/src/app/* ]]; then
      output_in_container="$output_path"
    elif [[ "$first_input" == /usr/src/app/* ]]; then
      output_in_container="/usr/src/app/$output_path"
    else
      output_in_container="/work/$output_path"
    fi

    mapped_args+=("$command_name" "$output_in_container" "${mapped_inputs[@]}")
    if [[ $# -gt 0 ]]; then
      mapped_args+=("$@")
    fi
    ;;

  merge-columns)
    if [[ $# -lt 3 ]]; then
      echo "error: merge-columns requires <base-file> <extra-file> <output-file>" >&2
      exit 2
    fi

    base_path="$1"
    extra_path="$2"
    output_path="$3"
    shift 3

    base_in_container="$(map_to_container_path "$base_path")"
    extra_in_container="$(map_to_container_path "$extra_path")"
    if [[ "$output_path" == /work/* || "$output_path" == /usr/src/app/* ]]; then
      output_in_container="$output_path"
    elif [[ "$base_in_container" == /usr/src/app/* ]]; then
      output_in_container="/usr/src/app/$output_path"
    else
      output_in_container="/work/$output_path"
    fi

    mapped_args+=("$command_name" "$base_in_container" "$extra_in_container" "$output_in_container")
    if [[ $# -gt 0 ]]; then
      mapped_args+=("$@")
    fi
    ;;

  reorder-columns)
    if [[ $# -lt 2 ]]; then
      echo "error: reorder-columns requires <input-file> <output-file>" >&2
      exit 2
    fi
    input_path="$1"
    output_path="$2"
    shift 2

    input_in_container="$(map_to_container_path "$input_path")"
    if [[ "$output_path" == /work/* || "$output_path" == /usr/src/app/* ]]; then
      output_in_container="$output_path"
    elif [[ "$input_in_container" == /usr/src/app/* ]]; then
      output_in_container="/usr/src/app/$output_path"
    else
      output_in_container="/work/$output_path"
    fi

    mapped_args+=("$command_name" "$input_in_container" "$output_in_container")
    if [[ $# -gt 0 ]]; then
      mapped_args+=("$@")
    fi
    ;;

  *)
    echo "error: unknown command '$command_name' (expected: filter-columns, filter-rows, replace-rows, merge-columns, reorder-columns, or compare-columns)" >&2
    exit 2
    ;;
esac

cmd=(docker compose --profile tools run)
if [[ -n "$compose_build_flag" ]]; then
  cmd+=("$compose_build_flag")
fi
if [[ -n "$compose_no_cache_flag" ]]; then
  cmd+=("$compose_no_cache_flag")
fi
cmd+=(--rm utilities)
cmd+=("${mapped_args[@]}")

exec "${cmd[@]}"
