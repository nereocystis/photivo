#ifndef FINALLY_H
#define FINALLY_H

#include <functional>

class finally
{
public:
    explicit finally(std::function<void()> f)
        : finally_clause(f)
    {
    }

    void reset()
    {
        finally_clause = std::function<void()>();
    }

    ~finally()
    {
        if (finally_clause)
        {
            finally_clause();
        }
    }

private:
    std::function<void()> finally_clause;

    finally(const finally &);
    finally &operator=(const finally &);
};

#endif /* FINALLY_H */
