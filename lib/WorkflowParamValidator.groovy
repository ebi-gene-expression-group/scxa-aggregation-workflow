import java.util.regex.Pattern

class WorkflowParamValidator {
    private static final Pattern INTEGER = Pattern.compile(/[0-9]+/)
    private static final Pattern PATH_VALUE = Pattern.compile(/[A-Za-z0-9._+,:=@%\/-]+/)
    private static final Set LEVELS = ['gene', 'transcript'] as Set
    private static final Set SCALING_MODES = ['no', 'scaledTPM', 'lengthScaledTPM', 'dtuScaledTPM'] as Set
    private static final Set BOOLEAN_STRINGS = ['TRUE', 'FALSE'] as Set

    static void validate(def params) {
        requirePath(params, 'resultsRoot')
        requirePath(params, 'quantDir')
        requireEnum(params, 'level', LEVELS)
        requireEnum(params, 'scaling', SCALING_MODES)
        requireInteger(params, 'chunkSize')
        requireNestedEnum(params.reference, 'params.reference', 'ignoreTxVersion', BOOLEAN_STRINGS)
    }

    private static void requirePath(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", PATH_VALUE)
    }

    private static void requireInteger(def params, String name) {
        requireValue(params, "params.${name}", name)
        assertPattern(params.get(name), "params.${name}", INTEGER)
    }

    private static void requireEnum(def params, String name, Set allowed) {
        requireValue(params, "params.${name}", name)
        def text = params.get(name).toString()
        if (!(text in allowed)) {
            throw new IllegalArgumentException("params.${name} must be one of ${allowed}; got '${text}'")
        }
    }

    private static void requireNestedEnum(def params, String scope, String name, Set allowed) {
        requireValue(params, "${scope}.${name}", name)
        def text = params.get(name).toString()
        if (!(text in allowed)) {
            throw new IllegalArgumentException("${scope}.${name} must be one of ${allowed}; got '${text}'")
        }
    }

    private static void requireValue(def params, String label, String key) {
        if (params == null || !has(params, key) || params.get(key) == null || params.get(key).toString() == '') {
            throw new IllegalArgumentException("Missing required workflow parameter ${label}")
        }
    }

    private static void assertPattern(value, String label, Pattern pattern) {
        def text = value == null ? '' : value.toString()
        if (!pattern.matcher(text).matches()) {
            throw new IllegalArgumentException("Invalid workflow parameter ${label}: '${text}'")
        }
    }

    private static boolean has(def params, String key) {
        params != null && params.containsKey(key)
    }
}
