
public class smFRETAnalysisException extends RuntimeException {

    // Default constructor
    public smFRETAnalysisException() {
        super();
    }

    // Constructor that accepts a custom error message
    public smFRETAnalysisException(String message) {
        super(message);
    }

    // Constructor that chains another exception (preserves root cause)
    public smFRETAnalysisException(String message, Throwable cause) {
        super(message, cause);
    }
}
